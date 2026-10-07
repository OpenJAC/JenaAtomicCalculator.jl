

# The EOL (optimized-level) path: the CSF-pair coefficient cache, the CI matrix and its diagonalization,
# the generalized occupation, and the two EOL solvers.

# Said once per session, not once per layer: a ladder prints it on every step otherwise, and the note is
# advice about how the run was STARTED rather than about anything the layer is doing.
const GBL_EOL_THREAD_NOTE_SHOWN = Ref(false)


"""
`struct  SelfConsistent.ScfVerdict`
    ... what the EOL rotation solver concluded about its own run, in numbers a caller can assert on rather than
        text a caller must grep.  It is ADVICE and never a gate: nothing in JAC reads it, no computation is
        blocked by it, and a surprising physical result must remain obtainable with every field looking wrong.

    + converged            ::Bool      ... advice, NOT a proof.  It means the run reached a stationary point TO
                                           MACHINE RESOLUTION FROM THIS START.  It does NOT mean the minimum:
                                           measured 10-Sep-2026, two legitimate starts on one system land 0.05 to
                                           0.57 mHa apart and BOTH report converged.  Where a milli-Hartree
                                           matters, run from more than one start and keep the lowest.
    + stopReason           ::String    ... which exit fired, verbatim: "converged", "energy below its own
                                           resolution", "no descent", "energy stagnated", "direction collapsed",
                                           "all subshells frozen", or "" for the iteration budget.
    + iterations           ::Int64     ... iterations actually performed.
    + energyStillAvailable ::Float64   ... Ha still reachable along the search direction, the Newton decrement
                                           dg^2/(2 <d,Hd>).  **This is the field to assert on**, because |grad|
                                           cannot be compared across systems and this can: measured on one ladder,
                                           4.50e-11 Ha on a finished layer against 5.81e-06 on an unfinished one
                                           where |grad| differed by only 48x.  A LOWER BOUND -- only one direction
                                           is measured -- and NaN where it was not measured or the curvature along
                                           that direction was not positive.
    + gradientNorm         ::Float64   ... |grad| at the exit, kept as a HINT only.  It is not scale-free, and it
                                           is partly supported on rotations the solver is forbidden to make: the
                                           raw gradient carries negative-energy-branch components whose curvature
                                           is -1.5e+04 Ha.  Do not threshold on it.
    + finalStep            ::Float64   ... the line-search step at the exit.  A flat energy means two different
                                           things and only this tells them apart: flat with a healthy step is
                                           convergence, flat with a collapsed step is a search that cannot move.
    + energy               ::Float64   ... the active-part energy at the exit, for reference.
    + descent              ::Float64   ... how far this solve brought the energy DOWN from its own first iterate,
        E(1) - E(final), in Hartree.  Within one RAS layer the frozen part is constant, so this is the layer's
        INCREMENT -- the quantity Rule 21 says to record and compare, rather than the total.
    + descentAtHalf        ::Float64   ... the same descent measured at HALF the iterations it actually took,
        E(1) - E(n/2).  The pair is what makes the doubling test free: a correlation layer never truly converges,
        so `maxIterations` is a cost dial and not a criterion, and the honest question is whether DOUBLING the
        effort still moves the increment.  Comparing these two answers it from one run, where it used to take two
        -- and more correctly, since both numbers come from ONE trajectory rather than from two computations whose
        earlier layers might have stopped in different places.
"""
struct  ScfVerdict
    converged              ::Bool
    stopReason             ::String
    iterations             ::Int64
    energyStillAvailable   ::Float64
    gradientNorm           ::Float64
    finalStep              ::Float64
    energy                 ::Float64
    descent                ::Float64
    descentAtHalf          ::Float64
end


"""
`Base.show(io::IO, verdict::SelfConsistent.ScfVerdict)`
    ... prints the verdict of the last EOL rotation solve in one block; nothing is returned.
"""
function Base.show(io::IO, verdict::SelfConsistent.ScfVerdict)
    println(io, "ScfVerdict:  $(verdict.converged ? "converged" : "NOT converged") " *
                "($(verdict.stopReason == "" ? "iteration budget" : verdict.stopReason)) " *
                "after $(verdict.iterations) iterations")
    println(io, "   energy still available = $(verdict.energyStillAvailable) Ha   (lower bound; the field to " *
                "assert on)")
    println(io, "   |grad| = $(verdict.gradientNorm) (a HINT, not scale-free)   final step = $(verdict.finalStep)" *
                "   energy = $(verdict.energy)")
    # THE DOUBLING TEST, FOR FREE.  Both numbers come from the one trajectory, so this says what a second
    # computation at half the budget would have said -- see ScfVerdict's own documentation.
    if  verdict.descent > 0.  &&  verdict.descentAtHalf > 0.
        drift = abs(verdict.descent - verdict.descentAtHalf) / verdict.descent
        println(io, "   descent = $(verdict.descent) Ha, and $(verdict.descentAtHalf) Ha at half the iterations: " *
                    "doubling the effort moved it by " * @sprintf("%.2f %%", 100*drift) *
                    (drift <= 0.10 ? "  -- stable" : "  -- NOT stable, the increment is still moving"))
    end
end


# The last verdict, set by every exit of solveOptimizedLevelFieldByRotation and read by
# SelfConsistent.lastScfVerdict().  A global is sound here precisely because this solver is single-threaded --
# see the note in the cost-estimate block -- and it is what keeps the change NON-BREAKING: the solver still
# returns its Multiplet, so no caller has to be touched to gain a verdict it can test.
const GBL_EOL_LAST_VERDICT = Ref{Union{Nothing, ScfVerdict}}(nothing)


"""
`SelfConsistent.lastScfVerdict()`
    ... returns what the EOL rotation solver concluded about its most recent run, so that a script can ASK whether
        a computation converged instead of capturing stdout and grepping it.  A `verdict::ScfVerdict` is returned,
        or `nothing` if that solver has not run in this session.

        MEASURED COST OF NOT HAVING THIS, and the reason the item existed: on 27-Sep-2026 a tuning scan suppressed
        the SCF output to keep a table readable, and an UNCONVERGED ground-state SCF then looked like a physics
        result -- the excitation energy came out non-monotonic and was one step from being written up as a
        physical gate. It was an artefact, and it was caught only by re-running with the log captured and the
        markers counted by hand.

        IT IS ADVICE AND NEVER A GATE. Read `energyStillAvailable` rather than `converged` wherever a number will
        do, and read `ScfVerdict`'s own documentation for what `converged` does and does not claim.
"""
function lastScfVerdict()
    return( GBL_EOL_LAST_VERDICT[] )
end

"""
`struct  SelfConsistent.PairCoefficientCache`
    ... holds the orbital-independent angular coefficients of every CSF pair of ONE symmetry block, in a form whose
        size is set by the number of DISTINCT labels and values rather than by the number of pair-entries.

        WHY IT IS NOT A `Dict` OF `Vector`s ANY MORE.  It was, until 14-Sep-2026, and that cost 1.55 MB per CSF
        against 0.181 for a plain Coulomb CI of the same ion -- about ten times more -- which is what put an open
        d shell out of reach whatever the machine: a Ti III layer of 25 085 CSF was killed at a predicted 72 GB.
        The CI matrix was never the problem; it IS built and diagonalised per symmetry block. The STORAGE was,
        because every block's cache is held for the whole SCF run (the EOL target levels span symmetries and each
        outer iteration revisits them) and because the old form paid a `Vector` header for every one of the n^2
        ordered pairs, including the empty ones.

        THE THREE REDUNDANCIES, MEASURED 13-Sep-2026 on the J = 2+ block of a Ti III 3p-opened space, 239 CSF:
            57 121 ordered pairs, of which 34 686 (61 %) carry NOTHING and still held an empty Vector;
            140 585 coefficient entries drawn from only 437 distinct (nu,a,b,c,d) labels   -> 321 x reuse;
            the same entries drawn from 15 560 distinct values                             ->   9 x reuse.
        So: FLAT CSR arrays instead of a Dict of Vectors, so an empty pair costs one Int32 and no header; the
        labels INTERNED in a table, an Int32 index per entry; the values likewise. Measured 11.4 x smaller on
        that block, and the gain GROWS with block size, which is where it is needed.

        NOT DONE, and deliberately: the further factor of two from hermiticity. `buildCIMatrixEOL` uses only
        r <= s, but `combineAngularCoefficientsEOL` uses every ordered pair, so dropping the lower triangle needs
        the relation between the (r,s) and (s,r) coefficient lists PROVED rather than assumed. It is the smallest
        of the available gains and the only one resting on an unverified claim.

    + idxCsf     ::Array{Int64,1}            ... indices, in the full basis, of the CSFs of this block.
    + subshells  ::Array{Subshell,1}         ... the basis subshell list, so a label index can be turned back
                                                 into a Subshell.
    + labels1p   ::Array{NTuple{3,Int64},1}  ... distinct (nu, a, b) as subshell INDICES.
    + labels2p   ::Array{NTuple{5,Int64},1}  ... distinct (nu, a, b, c, d) as subshell indices.
    + values     ::Array{Float64,1}          ... distinct coefficient values, shared by the one- and two-particle
                                                 parts.
    + ptr1p, lab1p, val1p                    ... CSR row pointer and the two index arrays, one-particle.
    + ptr2p, lab2p, val2p                    ... the same, two-particle.
"""
struct  PairCoefficientCache
    idxCsf              ::Array{Int64,1}
    subshells           ::Array{Subshell,1}
    labels1p            ::Array{NTuple{3,Int64},1}
    labels2p            ::Array{NTuple{5,Int64},1}
    values              ::Array{Float64,1}
    ptr1p               ::Array{Int32,1}
    lab1p               ::Array{Int32,1}
    val1p               ::Array{Int32,1}
    ptr2p               ::Array{Int32,1}
    lab2p               ::Array{Int32,1}
    val2p               ::Array{Int32,1}
end


"""
`SelfConsistent.coefficients1p(c::SelfConsistent.PairCoefficientCache, r::Int64, s::Int64)`
    ... rebuilds, lazily, the one-particle coefficients of the CSF pair (r,s) from the interned tables. A
        generator of `Coefficient1p` is returned, so nothing is allocated per pair and the caller's
        `for cf in ...` loop is unchanged.
"""
function coefficients1p(c::PairCoefficientCache, r::Int64, s::Int64)
    k = (s - 1) * length(c.idxCsf) + r
    return( ( Coefficient1p( c.labels1p[c.lab1p[i]][1], c.subshells[c.labels1p[c.lab1p[i]][2]],
                             c.subshells[c.labels1p[c.lab1p[i]][3]], c.values[c.val1p[i]] )
              for i = c.ptr1p[k]:(c.ptr1p[k+1] - 1) ) )
end


"""
`SelfConsistent.coefficients2p(c::SelfConsistent.PairCoefficientCache, r::Int64, s::Int64)`
    ... as `SelfConsistent.coefficients1p`, for the two-particle coefficients. A generator of `Coefficient2p`
        is returned.
"""
function coefficients2p(c::PairCoefficientCache, r::Int64, s::Int64)
    k = (s - 1) * length(c.idxCsf) + r
    return( ( Coefficient2p( c.labels2p[c.lab2p[i]][1], c.subshells[c.labels2p[c.lab2p[i]][2]],
                             c.subshells[c.labels2p[c.lab2p[i]][3]], c.subshells[c.labels2p[c.lab2p[i]][4]],
                             c.subshells[c.labels2p[c.lab2p[i]][5]], c.values[c.val2p[i]] )
              for i = c.ptr2p[k]:(c.ptr2p[k+1] - 1) ) )
end


"""
`SelfConsistent.cacheCsfPairCoefficientsEOL(sym::LevelSymmetry, basis::Basis)`
    ... computes, once, the orbital-independent angular coefficients of every CSF pair of the symmetry block
        `sym` and stores them interned. The coefficients depend only on the CSFs and the subshell list, never on
        the radial functions, so they survive every outer SCF+CI iteration unchanged -- which is the whole point
        of caching them. A `cache::SelfConsistent.PairCoefficientCache` is returned; see its docstring for the
        storage form and for what it replaced.
"""
function cacheCsfPairCoefficientsEOL(sym::LevelSymmetry, basis::Basis)
    idxCsf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J == sym.J   &&   basis.csfs[idx].parity == sym.parity    push!(idxCsf, idx)    end
    end
    n     = length(idxCsf)
    subIx = Dict( sh => i  for (i, sh) in enumerate(basis.subshells) )
    lab1D = Dict{NTuple{3,Int64},Int32}();    labels1p = NTuple{3,Int64}[]
    lab2D = Dict{NTuple{5,Int64},Int32}();    labels2p = NTuple{5,Int64}[]
    valD  = Dict{Float64,Int32}();            values   = Float64[]
    intern!(d, tab, key) = get!(d, key) do;   push!(tab, key);   Int32(length(tab))   end
    ptr1p = Int32[1];   lab1p = Int32[];   val1p = Int32[]
    ptr2p = Int32[1];   lab2p = Int32[];   val2p = Int32[]
    # THE COLUMN-MAJOR ORDER (s outer, r inner) matches the linear index the accessors form, k = (s-1)n + r.
    for  s = 1:n
        for  r = 1:n
            csfR = basis.csfs[idxCsf[r]];   csfS = basis.csfs[idxCsf[s]]
            for  cf  in  SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus),
                                                               csfR, csfS, basis.subshells)
                push!(lab1p, intern!(lab1D, labels1p, (cf.nu, subIx[cf.a], subIx[cf.b])))
                push!(val1p, intern!(valD,  values,   cf.T))
            end
            push!(ptr1p, Int32(length(lab1p) + 1))
            for  cf  in  SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, Basics.plus),
                                                         csfR, csfS, basis.subshells)
                push!(lab2p, intern!(lab2D, labels2p, (cf.nu, subIx[cf.a], subIx[cf.b], subIx[cf.c], subIx[cf.d])))
                push!(val2p, intern!(valD,  values,   cf.V))
            end
            push!(ptr2p, Int32(length(lab2p) + 1))
        end
    end

    return( PairCoefficientCache(idxCsf, basis.subshells, labels1p, labels2p, values,
                                 ptr1p, lab1p, val1p, ptr2p, lab2p, val2p) )
end


"""
`SelfConsistent.buildCIMatrixEOL(cache::SelfConsistent.PairCoefficientCache, orbitals::Dict{Subshell, Orbital},
                                 grid::Radial.Grid, potential::Radial.Potential)`
    ... (re-) builds the CI Hamiltonian matrix for one symmetry block from CACHED, orbital-independent
        angular coefficients (see cacheCsfPairCoefficientsEOL) and the CURRENT radial functions in
        orbitals -- algebraically identical to Hamiltonian.setupMatrixKinkAware's pure-Coulomb contribution
        (kink-aware InteractionStrength.XL_CoulombKinkAware, matching the rest of the average-level (ALField) line; Breit
        and QED are added only once, at the final Hamiltonian.performCIKinkAware call, exactly as for the AL
        scheme), but without repeating the (unchanged) angular-coefficient computation on every outer
        SCF+CI iteration. The trailing radial1pCache/radial2pCache arguments memoize each radial integral by
        its bare subshell labels (never by CSF-pair index), so a radial integral shared by many different
        CSF pairs -- e.g. the same "1s-1s" self-interaction appearing in every CSF's diagonal term -- is
        evaluated once per outer iteration rather than once per (r,s) occurrence; pass the SAME two Dicts
        into every block diagonalized within one outer iteration to also share the cache across blocks
        (found by profiling to be the dominant EOL cost for multi-CSF cases -- see
        project_eol_implementation.md). A  matrix::Array{Float64,2}  is returned.
"""
function buildCIMatrixEOL(cache::PairCoefficientCache, orbitals::Dict{Subshell, Orbital},
                          grid::Radial.Grid, potential::Radial.Potential,
                          radial1pCache::Dict{Tuple{Subshell,Subshell},Float64}          = Dict{Tuple{Subshell,Subshell},Float64}(),
                          radial2pCache::Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64} =
                                        Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}(),
                          vkCache::Dict{Tuple{Int64,Subshell,Subshell,Int64},Vector{Float64}} = Dict{Tuple{Int64,Subshell,Subshell,Int64},Vector{Float64}}())
    # radial1pCache/radial2pCache: keyed purely by subshell labels (never by CSF-pair index), so a radial
    # integral shared by MANY different CSF pairs -- e.g. the same "1s-1s" self-interaction appearing in
    # every CSF's diagonal term -- is evaluated once per outer SCF+CI iteration and reused, instead of once
    # per (r,s) occurrence. Pass the SAME cache Dicts in across every block diagonalized within one outer
    # iteration (they also depend only on the current orbitals, not on which block/CSF-pair references them);
    # a fresh empty cache per call (the default) is still correct, just without the cross-block reuse.
    #   vkCache IS A THIRD CACHE AT A DIFFERENT GRAIN, and it is the one that matters most here.  radial2pCache
    # deduplicates the finished integral R^nu(abcd) by its QUINTUPLE, which is the right key for a caller that
    # meets the same quadruple twice.  It cannot see the redundancy INSIDE those integrals: the expensive part of
    # R^nu(abcd) is the screened potential V_nu[b,d], a function of only TWO of the four orbitals, so every
    # distinct (a,c) sharing one (nu,b,d) rebuilt it from scratch.  Measured on C-like uranium, 20 083
    # coefficients carried 20 083 distinct quintuples and only 781 distinct (nu,b,d) triples -- a redundancy of
    # 25.7x that a quintuple key cannot reach.  Counted again 05-Oct-2026 on a Ca+ ladder (priority item 47),
    # `buildScreenedPotential` was entered 32 737 times for 125 distinct (nu,b,d), and 42.7 % of ALL the
    # allocation attributed to that routine arrived through THIS call site.
    #   WHAT IT IS WORTH, measured A/B on one build and closing priority item 47 on 05-Oct-2026.  C-like uranium,
    # 1s^2 2s^2 2p^2 with SD into {3s..4f}, 864 CSFs over 16 subshells:
    #       builds       57 423 -> 17 635        69.3 % fewer
    #       allocation    23.27 -> 16.71 GB      28.2 % less
    #       wall clock   1.72x / 1.96x / 1.87x   three pairs with the ORDER ALTERNATED, so ~1.9x
    #       energy       -14350.388262411614051 both ways, dE = 0.000e+00 -- BIT-IDENTICAL
    # The gain grows with BOTH the CSF count and the subshell count: 30.7 % fewer builds at 19 CSFs, 29.1 % at
    # 149, 55.2 % at 870, 69.3 % at 864-over-16-subshells.  A single contended pair had read 2.4x and would have
    # overstated it by a quarter, which is why the repeats are alternated.
    #   AND THIS IS THE COMPANION OF THE 4.04x OF 03-Sep-2026, which did exactly this for the EOL FIELD path and
    # left the CI MATRIX path untouched for a month -- hidden behind `radial2pCache`, which keys on the finished
    # integral's quintuple and therefore cannot see redundancy inside it.  If a third path is ever added, ask
    # which grain its cache keys on before assuming it is covered.
    #   The key carries the extent as well as (nu,b,d), and it is the caller's: it must not outlive the orbitals
    # it was built from, so it is created beside radial1p/radial2p and dies with them at the end of the outer
    # iteration.
    # Hermitian-symmetry shortcut (28-Jul-2026): only the UPPER triangle (r<=s) is computed -- exact, not
    # just safe, since this matrix feeds diagonalizeBlockEOL -> Basics.diagonalize(MatrixWithLinearAlgebra(),
    # ...), whose Symmetric(matrix) wrapper (default uplo=:U) already discards the lower triangle. See
    # Hamiltonian.setupMatrix's identical note for the confirming test.
    # JAC_NO_VKCACHE restores the pre-`cebc08f3` behaviour -- the UNCACHED method, one screened potential per
    # quadruple -- so that the saving can be measured A/B on ONE build rather than against a number remembered
    # from another tree.  It exists for that measurement and for re-taking it on a new system; it is not a tuning
    # knob and nothing reads it by default.  Read ONCE here rather than per coefficient: an ENV lookup is a Dict
    # probe and does not belong in a numerical inner loop.
    noVk = haskey(ENV, "JAC_NO_VKCACHE")
    n = length(cache.idxCsf);   matrix = zeros(Float64, n, n)
    for  r = 1:n
        for  s = r:n
            me = 0.
            for  cf in coefficients1p(cache, r, s)
                I_ab = get!(radial1pCache, (cf.a,cf.b)) do
                    RadialIntegrals.GrantIab(orbitals[cf.a], orbitals[cf.b], grid, potential)
                end
                me = me + cf.T * I_ab
            end
            for  cf in coefficients2p(cache, r, s)
                R_abcd = get!(radial2pCache, (cf.nu,cf.a,cf.b,cf.c,cf.d)) do
                    if  noVk
                        InteractionStrength.XL_CoulombKinkAware(cf.nu, orbitals[cf.a], orbitals[cf.b],
                                                                orbitals[cf.c], orbitals[cf.d], grid)
                    else
                        InteractionStrength.XL_CoulombKinkAware(cf.nu, orbitals[cf.a], orbitals[cf.b],
                                                                orbitals[cf.c], orbitals[cf.d], grid, vkCache)
                    end
                end
                me = me + cf.V * R_abcd
            end
            matrix[r,s] = me
        end
    end

    return( matrix )
end


"""
`SelfConsistent.diagonalizeBlockEOL(sym::LevelSymmetry, idxCsf::Array{Int64,1}, matrix::Array{Float64,2}, basis::Basis)`
    ... diagonalizes the (already-built) CI matrix of one symmetry block and reassigns the eigenvectors to
        Level instances w.r.t. the full basis (zero-padded outside this block), exactly as
        Hamiltonian.performCI/performCIKinkAware do internally per symmetry block -- factored out here so it
        can be called every EOL outer iteration without their Multiplet-merge/tabulate overhead.
        An  Array{Level,1}  is returned.
"""
function diagonalizeBlockEOL(sym::LevelSymmetry, idxCsf::Array{Int64,1}, matrix::Array{Float64,2}, basis::Basis)
    eigen  = Basics.diagonalize(MatrixWithLinearAlgebra(), matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        vector = zeros( length(basis.csfs) )
        for  (r, idx)  in  enumerate(idxCsf)    vector[idx] = eigen.vectors[ev][r]    end
        newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector )
        push!( levels, newlevel)
    end

    return( levels )
end


"""
`SelfConsistent.selectTargetLevelsEOL(mp::Multiplet, levelSelectionCI::LevelSelection)`
    ... determines the target level(s) for the EOL functional from an (energy-sorted) multiplet mp,
        using levelSelectionCI EXCLUSIVELY.  If configurations is given, the target set is selected by
        reference weight (see below) and, if symmetries is given as well, restricted to those J^P --
        the one-orbital-set-per-J route for RAS layers.  Otherwise either indices or symmetries may be
        given, never both.
        If symmetries is given, the target set is the LOWEST level of each listed symmetry (the classic
        EOL use case, e.g. the lowest of J=1/2^+ together with the lowest of J=3/2^+). If indices is
        given, the target set is those exact levels by their (global, energy-sorted) index. If
        levelSelectionCI is inactive or both arrays are empty, the target set defaults to the single
        lowest level overall -- a genuine OL (one-level) computation.
        An  Array{Level,1}  is returned.
"""
function referenceCsfIndices(basis::Basis, refConfigs::Array{Configuration,1})
    # A CSF belongs to a reference configuration when its NON-RELATIVISTIC occupations agree with that
    # configuration's shell by shell; the relativistic subshells of one shell are summed, so that a
    # configuration named as 5f^1 matches CSFs carrying 5f_5/2^1 or 5f_7/2^1 alike.
    refOccs = [ Dict{Shell,Int64}( sh => n  for (sh, n) in conf.shells  if n > 0 )  for conf in refConfigs ]
    indices = Int64[]
    for  (ic, csf)  in  enumerate(basis.csfs)
        occ = Dict{Shell,Int64}()
        for  (i, subsh)  in  enumerate(basis.subshells)
            if  csf.occupation[i] > 0
                sh = Shell(subsh.n, Basics.subshell_l(subsh));    occ[sh] = get(occ, sh, 0) + csf.occupation[i]
            end
        end
        if  any(r -> r == occ, refOccs)    push!(indices, ic)    end
    end
    return( indices )
end


function selectTargetLevelsEOL(mp::Multiplet, levelSelectionCI::LevelSelection)
    if  levelSelectionCI.active  &&  !isempty(levelSelectionCI.configurations)
        # SELECTION BY REFERENCE WEIGHT.  Selecting the EOL target set by index is not merely imprecise, it is
        # UNSTABLE: the indices refer to the energy-sorted multiplet, so the moment a correlation configuration
        # sinks below the reference -- which is what a doubly-excited layer routinely does before its orbitals
        # have settled -- indices 1..n stop pointing at the reference levels and the field begins optimizing the
        # intruders.  The functional then changes identity between iterations and the optimizer has no fixed
        # minimum to find; measured on Cf^17+ with an SD layer into {7s,7p}, the gradient plateaued at 0.69 while
        # the step collapsed to 5e-8, and the returned "ground state" was a 7p^2 level of the wrong J.
        # Weighting on the reference CSFs is immune to that, because it asks what a level IS and not where it
        # sits.  The premise is the ordinary one for this kind of work: the reference configurations are chosen
        # so that the level ORDER is already right and correlation only improves the energies.
        refIdx = SelfConsistent.referenceCsfIndices(mp.levels[1].basis, levelSelectionCI.configurations)
        if  isempty(refIdx)     error("SelfConsistent.selectTargetLevelsEOL(): none of the $(length(mp.levels[1].basis.csfs)) " *
                                      "CSFs in this basis belongs to the reference configurations " *
                                      "$(levelSelectionCI.configurations).")   end
        targetLevels = Level[];    best = 0.0
        for  level  in  mp.levels                              # energy-sorted
            w = sum( level.mc[r]^2  for r in refIdx );    best = max(best, w)
            if  w >= 0.5    push!(targetLevels, level)    end
        end
        # A THRESHOLD OF 0.5 IS THE STATEMENT "this level IS a reference level", the same rule the Os^16+ and
        # Cf application reports had to apply by hand afterwards.  Failing it is not something to paper over
        # with a fallback to the lowest levels: it means the reference space does not describe this spectrum,
        # and continuing would optimize on states nobody chose.
        if  isempty(targetLevels)   error("SelfConsistent.selectTargetLevelsEOL(): no level carries a weight of 0.5 " *
                    "or more on the reference configurations $(levelSelectionCI.configurations); the largest is " *
                    "$best.  The reference space does not describe these levels, so no EOL target set can be formed.")
        end
        # A SYMMETRY LIST MAY BE GIVEN IN ADDITION, and then restricts the reference levels to those J^P.  This is the
        # ONE-ORBITAL-SET-PER-J route, and it is the recommended way to run the RAS layers on top of an AL reference:
        # the spectroscopic orbitals are optimized once, by AL, over all fine-structure levels of the reference
        # configurations; each added correlation layer is then optimized SEPARATELY for each J, giving one set of
        # correlation orbitals per symmetry.  Statistically (2J+1)-weighting several J's into a single set of
        # correlation orbitals asks those orbitals for a compromise no single J wants, and the price is paid in the
        # ENERGY: measured on Cf^17+ with an SD layer into {7s,7p}, the J=5/2 reference level comes out at
        # -32975.77417097 Ha when four levels share one orbital set and at -32975.77423934 Ha when it is optimized
        # alone -- 6.8e-5 Ha, about 15 cm^-1, and the clock transition moves from 8834 to 8849 cm^-1 with it.  Each
        # per-J energy is a proper variational upper bound for ITS OWN level, which the weighted compromise is not.
        # (An earlier version of this note claimed the compromise also collapses the line search, on the strength of
        # the step falling to 1e-5 in the four-level run.  That was a misreading: the energy of that run had already
        # converged to 1e-11 by iteration 11, so the shrinking step is an optimizer standing still at a flat minimum,
        # which is normal terminal behaviour and not a trust-radius pathology.)
        # Note this is a RESTRICTION of the reference-weight set and not the plain-symmetries branch below: the level
        # is chosen by what it IS (its weight on the reference configurations) and only then filtered by its J^P, so
        # an intruder of the right symmetry sinking below it cannot capture the target set.
        if  !isempty(levelSelectionCI.symmetries)
            selected = [ lv  for lv in targetLevels  if  LevelSymmetry(lv.J, lv.parity) in levelSelectionCI.symmetries ]
            if  isempty(selected)   error("SelfConsistent.selectTargetLevelsEOL(): of the $(length(targetLevels)) " *
                        "reference levels found, symmetries $(unique([LevelSymmetry(lv.J, lv.parity) for lv in targetLevels])), " *
                        "none carries one of the requested symmetries $(levelSelectionCI.symmetries).")
            end
            targetLevels = selected
        end
        return( targetLevels )
    elseif  !levelSelectionCI.active  ||  ( isempty(levelSelectionCI.indices) && isempty(levelSelectionCI.symmetries) )
        return( [ mp.levels[1] ] )
    elseif  !isempty(levelSelectionCI.indices)  &&  !isempty(levelSelectionCI.symmetries)
        error("stop a; levelSelectionCI must specify EITHER indices OR symmetries for the EOL scheme, not both.")
    elseif  !isempty(levelSelectionCI.symmetries)
        targetLevels = Level[]
        for  sym  in  levelSelectionCI.symmetries
            for  level  in  mp.levels                                   # mp.levels is energy-sorted; first match = lowest
                if  LevelSymmetry(level.J, level.parity) == sym    push!(targetLevels, level);   break    end
            end
        end
        return( targetLevels )
    else
        return( [ mp.levels[i]  for i in levelSelectionCI.indices ] )
    end
end


"""
`SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, targetLevels::Array{Level,1}, basis::Basis)`
    ... computes the EOL generalized occupation number per subshell,
        q(nlj) = Σᵣ d²_r · q_r(nlj), with d²_r = Σᵢ WT_i · (c_r⁽ⁱ⁾)² the DIAGONAL (r=r) case of the same
        statistical-(2J+1)-weighted generalized weight used in combineAngularCoefficientsEOL, and q_r(nlj)
        the plain occupation of nlj in CSF r (basis.csfs[r].occupation). This REPLACES
        Basics.extractMeanOccupation(basis) -- which averages FLATLY over every CSF in the whole basis,
        appropriate only for AL's single-average-CSF philosophy -- with the occupation actually implied by
        the current target level(s)' own CI mixing, recomputed every outer iteration since the mixing
        coefficients change. Without this, computeFockMatrix's (1.0/occ) two-electron scaling uses a
        basis-wide average that can be wildly wrong for a multi-configuration target level (e.g. an
        essentially-pure 2s² level in a 2s²/2p² basis would otherwise see occ(2s) diluted by the unrelated
        2p² CSFs, inflating its two-electron potential many-fold). A  Dict{Subshell,Float64}  is returned.
"""
function computeGeneralizedOccupationEOL(blockCaches, targetLevels::Array{Level,1}, basis::Basis)
    twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
    sumWeights  = sum( twiceJp1(level.J)  for level in targetLevels )
    weights     = [ twiceJp1(level.J) / sumWeights  for level in targetLevels ]

    occs = Dict{Subshell, Float64}();   for  sh in basis.subshells   occs[sh] = 0.   end
    for  (_, cache)  in  blockCaches
        idxCsf = cache.idxCsf
        for  r  in  idxCsf
            drr = 0.
            for  (i, level)  in  enumerate(targetLevels)    drr = drr + weights[i] * level.mc[r]^2    end
            if  drr == 0.    continue    end
            for  (is, sh)  in  enumerate(basis.subshells)   occs[sh] = occs[sh] + drr * basis.csfs[r].occupation[is]   end
        end
    end

    return( occs )
end


"""
`SelfConsistent.computeOrbitalEnergiesEOL(subshells::Array{Subshell,1}, bVectors::Dict{Subshell, Vector{Float64}},
                                          coeffs2p::Array{Coefficient2p,1}, genOcc::Dict{Subshell,Float64},
                                          primitives::Bsplines.Primitives, grid::Radial.Grid, nucPot::Radial.Potential,
                                          storage::Dict{String,Array{Float64,2}}, matrixB::Array{Float64,2};
                                          coeffs2pUnscaled::Array{Coefficient2p,1}=Coefficient2p[])`
    ... gives every subshell of a converged EOL orbital set a DEFINED energy, the diagonal Lagrange multiplier

            eps_a  =  <a|F_a|a>  =  b_a^T F_a b_a / b_a^T B b_a ,

        with F_a the same orbital-specific Fock operator `SelfConsistent.computeFockMatrix` builds to refine that
        subshell. This is the exact analogue of what the AL field already stores: there the orbital IS an eigenvector
        of F_a, so its stored energy `wc.values[ni]` and this Rayleigh quotient are the same number. A
        rotation-optimized orbital is NOT such an eigenvector, and the quotient is then the honest generalisation --
        the multiplier that the orthonormality constraint carries -- rather than an eigenvalue it does not have.

        WHAT THIS IS NOT: the one-particle expectation <a|h_D|a>. That omits the electron-electron interaction
        entirely and comes out 1.7x to 12x too DEEP (measured 12-Sep-2026, `tools/probe-eolOrbitalEnergy.jl`, Be-like
        and Ne-like: 2s of Ne-like at -11.13 against the AL -1.94), while looking perfectly well-behaved. It is the
        plausible wrong answer and must not be substituted here.

        A subshell carrying ZERO generalized occupation keeps 0.0, and that is correct rather than a fallback: the
        EOL functional does not depend on such an orbital at all, so no mean field and no multiplier is defined for
        it. A `Dict{Subshell, Float64}` is returned.
"""
function computeOrbitalEnergiesEOL(subshells::Array{Subshell,1}, bVectors::Dict{Subshell, Vector{Float64}},
                                    coeffs2p::Array{Coefficient2p,1}, genOcc::Dict{Subshell,Float64},
                                    primitives::Bsplines.Primitives, grid::Radial.Grid, nucPot::Radial.Potential,
                                    storage::Dict{String,Array{Float64,2}}, matrixB::Array{Float64,2};
                                    coeffs2pUnscaled::Array{Coefficient2p,1}=Coefficient2p[])
    energies    = Dict{Subshell, Float64}()
    neededRanks = unique( [ cf.nu  for cf in vcat(coeffs2p, coeffs2pUnscaled) ] )
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
        tensorCaches[L] = (cacheLL, cacheLS, cacheSS)
    end
    directKernels   = Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}}()
    exchangeKernels = Dict{Tuple{Int64,Subshell},Array{Float64,2}}()

    for  sh  in  subshells
        occ = get(genOcc, sh, 0.0)
        if  abs(occ) < 1.0e-12    energies[sh] = 0.0;    continue    end
        matrix = SelfConsistent.computeFockMatrix(sh, coeffs2p, bVectors, primitives, nucPot, storage, occ,
                                                  tensorCaches; coeffs2pUnscaled=coeffs2pUnscaled,
                                                  directKernels=directKernels, exchangeKernels=exchangeKernels)
        b            = bVectors[sh]
        energies[sh] = (transpose(b) * matrix * b) / (transpose(b) * matrixB * b)
    end

    return( energies )
end


"""
`SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels::Array{Level,1})`
    ... generalizes SelfConsistent.computeAngularCoefficients (AL's single-CSF-average analog: loop CSFs,
        weight 1/ncsf) to CSF PAIRS, weighted by the EOL generalized weight
        d²_rs = Σᵢ WT_i · c_r⁽ⁱ⁾ · c_s⁽ⁱ⁾, with WT_i the normalized statistical (2Jᵢ+1) weight of each
        target level i in targetLevels (its own mixing vector c⁽ⁱ⁾ = level.mc) -- never a user-supplied
        weight. Reuses the cached, orbital-independent per-CSF-pair coefficients from
        cacheCsfPairCoefficientsEOL for every relevant symmetry block in blockCaches; a level's mc vector
        is zero outside its own block, so multiple symmetry blocks combine correctly without
        special-casing. The dedup/condensation logic is identical to computeAngularCoefficients.
        A Tuple  (coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1})  is returned.
"""
function combineAngularCoefficientsEOL(blockCaches, targetLevels::Array{Level,1}; pairs::Symbol=:all)
    twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
    sumWeights  = sum( twiceJp1(level.J)  for level in targetLevels )
    weights     = [ twiceJp1(level.J) / sumWeights  for level in targetLevels ]

    # ACCUMULATED DIRECTLY INTO THE CONDENSED SET, never through a flat vector of every pair's coefficients.
    # It was built that way until 14-Sep-2026: the temporary held every coefficient of every contributing CSF
    # pair and was then condensed by a NESTED DOUBLE LOOP over that vector -- O(N^2) in the entry count, not in
    # the CSF count.  Accumulating directly removes both.  Measured on a 7 062-CSF Ti III space: 208.8 s ->
    # 130.1 s, a factor 1.6, with the energies 1 ulp apart.
    #
    # THE RESULT WAS NEVER LARGE: it is bounded by the number of distinct (nu,a,b,c,d) labels, i.e. by the
    # SUBSHELL count, while the temporary was bounded by nothing.
    #
    # TWO CLAIMS THAT LED HERE WERE LATER CHECKED AND ARE FALSE.  They are recorded because both are the kind a
    # reader would otherwise re-derive, and one of them stood in this comment until 30-Sep-2026.
    #   * NOT A BOXING PROBLEM.  The vector was said to hold separately boxed objects because `Coefficient2p[]`
    #     has a non-concrete element type.  It does not: inside this module `Coefficient1p` and `Coefficient2p`
    #     are `const` aliases to the CONCRETE types (module-SelfConsistent.jl, just below the Rule 18 note), so
    #     those vectors have always been concretely typed.  The "105-178 bytes per entry" that started it was an
    #     average over one-particle (48 B) and two-particle (80 B) entries plus array overhead, divided by the
    #     combined count; a benchmark that appeared to confirm it had used `SpinAngular.Coefficient1p`, the
    #     UnionAll, rather than the alias.
    #   * AND IT DOES NOT SET THE PEAK.  A VmHWM attribution gave this routine 0.584 GB of an EOL run's
    #     high-water mark, which is what made it look like the dominant term.  Removing the temporary saved
    #     158 MB, not 584, because in a full solve the peak is set later in the iteration.  The memory of a
    #     large EOL run is ALLOCATION CHURN -- measured 133.5 GB allocated to hold ~2 GB live, 12 % of the wall
    #     clock in collection -- so it is a rate, not a structure, and no single routine owns it.
    #
    # FIRST-ENCOUNTER ORDER IS PRESERVED DELIBERATELY, and it is what makes the change verifiable: the old code
    # summed each label's contributions in the order it met them and emitted the labels in that same order, so
    # keeping both leaves every downstream floating-point sum in the order it had.  Emitting the table in Dict
    # order instead would re-order those sums and move the last bits of the energies.
    idx1 = Dict{NTuple{3,Any},Int64}();    keys1 = NTuple{3,Any}[];    vals1 = Float64[]
    idx2 = Dict{NTuple{5,Any},Int64}();    keys2 = NTuple{5,Any}[];    vals2 = Float64[]
    for  (_, cache)  in  blockCaches
        idxCsf = cache.idxCsf;    n = length(idxCsf)
        for  r = 1:n
            for  s = 1:n
                # pairs = :all (default, unchanged) | :diagonal (r == s only) | :offdiagonal (r != s only).
                if      pairs == :diagonal      &&  r != s     continue
                elseif  pairs == :offdiagonal   &&  r == s     continue
                end
                drs = 0.
                for  (i, level)  in  enumerate(targetLevels)    drs = drs + weights[i] * level.mc[idxCsf[r]] * level.mc[idxCsf[s]]    end
                if  drs == 0.    continue    end
                for  cf in coefficients1p(cache, r, s)
                    key = (cf.nu, cf.a, cf.b)
                    ix  = get(idx1, key, 0)
                    if  ix == 0   push!(keys1, key);   push!(vals1, cf.T * drs);   idx1[key] = length(keys1)
                    else          vals1[ix] = vals1[ix] + cf.T * drs
                    end
                end
                for  cf in coefficients2p(cache, r, s)
                    key = (cf.nu, cf.a, cf.b, cf.c, cf.d)
                    ix  = get(idx2, key, 0)
                    if  ix == 0   push!(keys2, key);   push!(vals2, cf.V * drs);   idx2[key] = length(keys2)
                    else          vals2[ix] = vals2[ix] + cf.V * drs
                    end
                end
            end
        end
    end

    coeffs1px = Coefficient1p[];     coeffs2px = Coefficient2p[]
    for  (i, key)  in  enumerate(keys1)    push!(coeffs1px, Coefficient1p(key[1], key[2], key[3], vals1[i]))    end
    for  (i, key)  in  enumerate(keys2)
        push!(coeffs2px, Coefficient2p(key[1], key[2], key[3], key[4], key[5], vals2[i]))
    end

    return( (coeffs1px, coeffs2px) )
end


"""
`SelfConsistent.energyFromBVectors(bVectors::Dict{Subshell, Vector{Float64}},
        coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, grid::Radial.Grid,
        nucPot::Radial.Potential)`  
    ... the EOL energy as a plain scalar function of the orbital B-spline coefficient vectors, with the
        angular coefficients held fixed. This is exactly the functional solveOptimizedLevelField reports,
        just expressed in the variables that are actually varied, so that it can be differentiated. Note
        that it includes coeffs1p -- which the Fock matrix of the present scheme never receives, and which
        is the inconsistency the rotation-based path exists to remove. A value::Float64 is returned.
"""
function energyFromBVectors(bVectors::Dict{Subshell, Vector{Float64}},
                                   coeffs1p::Array{Coefficient1p,1},
                                   coeffs2p::Array{Coefficient2p,1},
                                   subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                   grid::Radial.Grid, nucPot::Radial.Potential)
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives; canonicalize=false)
    end
    return( SelfConsistent.computeFunctional(coeffs1p, coeffs2p, orbitals, grid, nucPot) )
end


"""
`SelfConsistent.energyFromBVectorsSplit(bVectors, coeffs1p, coeffs2p, subshells, primitives, grid, nucPot,
                                        isFrozen, frozenRk)`  
    ... as energyFromBVectors, but returns the pair (eFrozen, eActive) of SelfConsistent.computeFunctionalSplit,
        which is what the line search compares and what carries the persistent radial memo.
"""
function energyFromBVectorsSplit(bVectors::Dict{Subshell, Vector{Float64}},
                                 coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
                                 subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                 grid::Radial.Grid, nucPot::Radial.Potential,
                                 isFrozen::Function, frozenRk::Dict{NTuple{5,Any}, Float64})
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives; canonicalize=false)
    end
    return( SelfConsistent.computeFunctionalSplit(coeffs1p, coeffs2p, orbitals, grid, nucPot, isFrozen, frozenRk) )
end



"""
`SelfConsistent.positiveBranchSpectrum(subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
        nucPot::Radial.Potential, matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}})`  
    ... the positive-energy eigenvectors of the ONE-ELECTRON Dirac matrix of every kappa present, together
        with their eigenvalues. It depends only on the nuclear potential and the B-spline basis, NOT on the
        orbitals, so it is constant for a whole run and must be built once and passed in -- both
        virtualDirections and projectOntoPositiveBranch used to rebuild it on every call,
        i.e. twice per iteration. Measured on Be-like C^2+: an iteration cost 1.77 s of which the line
        search was only 0.15 s, and this was the bulk of the remainder.
        A Dict{Int64, Tuple{Array{Vector{Float64},1}, Vector{Float64}}} is returned, keyed by kappa.
"""
function positiveBranchSpectrum(subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       nucPot::Radial.Potential, matrixB::Array{Float64,2},
                                       storage::Dict{String,Array{Float64,2}})
    spectrum = Dict{Int64, Tuple{Array{Vector{Float64},1}, Vector{Float64}}}()
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        oneEl = Bsplines.setupLocalMatrix(kappa, primitives, nucPot, storage)
        wc    = Bsplines.diagonalizeLocalMatrix(kappa, oneEl, matrixB, primitives)
        mm    = Bsplines.findPositiveBranchStart(wc.values)
        spectrum[kappa] = ( [ wc.vectors[i]  for i = mm:length(wc.values) ],
                            [ wc.values[i]   for i = mm:length(wc.values) ] )
    end
    return( spectrum )
end


"""
`SelfConsistent.virtualDirections(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
        primitives::Bsplines.Primitives, nucPot::Radial.Potential, matrixB::Array{Float64,2},
        storage::Dict{String,Array{Float64,2}}; nVirtual::Int64=20)`  
    ... builds, for each subshell, an S-orthonormal set of ALLOWED rotation directions: positive-branch
        eigenvectors of the one-electron Dirac matrix of that kappa, projected free of every occupied
        orbital of the same kappa. Rotations among the occupied orbitals themselves are deliberately
        excluded -- they leave the CSF space invariant, are redundant, and would make the Hessian singular.
        A Dict{Subshell, Array{Vector{Float64},1}} is returned.
"""
function virtualDirections(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
                                  primitives::Bsplines.Primitives, nucPot::Radial.Potential,
                                  matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}};
                                  nVirtual::Int64=20, spectrum=nothing)
    posSpec = isnothing(spectrum) ?
              SelfConsistent.positiveBranchSpectrum(subshells, primitives, nucPot, matrixB, storage) : spectrum
    virtuals = Dict{Subshell, Array{Vector{Float64},1}}()
    for  sh  in  subshells
        # the occupied orbitals of this kappa, S-orthonormalized so that the projection below is exact
        occSame = Vector{Float64}[]
        for  s2  in  subshells
            if  s2.kappa != sh.kappa    continue    end
            v = copy(bVectors[s2])
            for  u in occSame    v = v - (transpose(u) * matrixB * v) * u    end
            nrm = sqrt( abs(transpose(v) * matrixB * v) )
            if  nrm > 1.0e-10    push!(occSame, v / nrm)    end
        end
        # a one-electron reference spectrum for this kappa; orbital-independent, hence a fixed frame
        (posVecs, _) = posSpec[sh.kappa]
        dirs  = Vector{Float64}[];    nrms = Float64[]
        for  i = 1:length(posVecs)
            v = copy(posVecs[i])
            for  u in occSame    v = v - (transpose(u) * matrixB * v) * u    end
            for  u in dirs       v = v - (transpose(u) * matrixB * v) * u    end
            nrm = sqrt( abs(transpose(v) * matrixB * v) )
            # THE THRESHOLD IS 0.1 AND NOT 1e-8 -- 05-Sep-2026.  `nrm` is what SURVIVES the projection out of the
            # occupied space, so nrm^2 is the fraction of this reference vector that is genuinely new.  At 1e-8 a
            # vector 99.999999 % inside the occupied space was accepted and then DIVIDED BY that residue,
            # amplifying round-off and leftover occupied character by up to 1e8.  Measured on Be-like U, the
            # accepted residues per kappa fall into two groups with nothing between them:
            #     s_1/2 : 3.9e-03  1.3e-02 | 4.3e-01 7.2e-01 7.3e-01 ... 8.2e-01
            #     p_1/2 : 3.3e-02          | 5.7e-01 8.2e-01 ...
            # The first group is amplified 77x to 256x and is where the direction lying BELOW the occupied
            # orbital came from -- the one that carried 100 % of the step's weight through the +0.05 floor.
            # 0.1 sits in the gap and states a criterion: at least 1 % of the vector must lie outside the
            # occupied space for it to count as a direction.  Fifteen or so of sixteen survive.
            if  nrm > 0.1    push!(dirs, v / nrm);   push!(nrms, nrm)    end
            if  length(dirs) >= nVirtual    break    end
        end
        if  haskey(ENV, "JAC_EOL_VIRTCHECK")
            println(">> [EOL-NRM] $sh  residual norms of the accepted directions: " *
                    join([@sprintf("%.1e", x) for x in nrms], " "));   flush(stdout)
        end

        # AND NOW MAKE THEM ORBITALS AGAIN -- priority item 5, 04-Sep-2026.
        # Up to here `dirs` spans the right space but its members are MIXTURES: each reference eigenvector is
        # projected free of the occupied orbitals and then Gram-Schmidt'ed against the ones already accepted, and
        # both steps mix.  The caller divides the gradient by (epsV - epsA) with epsV = phi' h1 phi, which is the
        # standard orbital-energy denominator and is only meaningful when phi is an EIGENSTATE; for a mixture it
        # is a Rayleigh quotient of nothing in particular.  MEASURED CONSEQUENCE (Be-like Z = 92, the sixteen
        # 2s_1/2 directions): six lay BELOW the occupied 2s at -1256 Ha, the lowest beneath even the 1s, the
        # floor then replaced their large NEGATIVE denominator by the smallest positive one available (+0.05),
        # amplifying them and reversing their sign, and those six carried 100.0 % of the step's weight while
        # 2p_1/2 and 2p_3/2 carried none.
        # THE REPAIR IS TO DIAGONALIZE h1 INSIDE THE SPAN, which costs one nVirtual x nVirtual symmetric
        # eigenproblem per kappa and turns the mixtures back into eigenstates of the one-particle operator
        # RESTRICTED TO THE COMPLEMENT of the occupied space.  Each epsV the caller then forms is a genuine
        # eigenvalue rather than a Rayleigh quotient, which is exactly what its formula assumes.
        # TWO PROPERTIES ARE PRESERVED, and the callers rely on both: the set stays B-ORTHONORMAL, because an
        # orthogonal transformation of a B-orthonormal set is B-orthonormal; and it stays inside the complement,
        # because a linear combination of vectors orthogonal to the occupied space still is.
        # h1 IS THE ONE-PARTICLE OPERATOR, NOT THE FULL MEAN FIELD, deliberately: the caller's denominator uses
        # h1 on BOTH sides (epsA = b' h1 b), the rotation route never forms a Fock matrix, and a denominator is
        # self-consistent or it is nothing.  Eigenstates of the true mean field would be a larger change and are
        # not needed to make the formula mean what it says.
        if  length(dirs) > 1
            h1  = Bsplines.setupLocalMatrix(sh.kappa, primitives, nucPot, storage)
            nd  = length(dirs)
            hm  = zeros(nd, nd)
            for  i = 1:nd
                hv = h1 * dirs[i]
                for  j = i:nd    hm[i,j] = transpose(dirs[j]) * hv;    hm[j,i] = hm[i,j]    end
            end
            wa   = LinearAlgebra.eigen( LinearAlgebra.Symmetric(hm) )
            dirs = [ sum( wa.vectors[k,m] * dirs[k]  for k = 1:nd )   for m = 1:nd ]
        end
        virtuals[sh] = dirs
    end
    return( virtuals )
end


"""
`SelfConsistent.expandBVector(vec::Vector{Float64}, primitives::Bsplines.Primitives)`  
    ... expands a B-spline coefficient vector into its large and small radial components on the grid,
        (P, Q), exactly as the matrix contraction implies -- P(r) = sum_k vec[k] B_k^L(r) and likewise for Q.
        Deliberately NOT via Bsplines.generateOrbitalFromVector, which truncates at mtp and cleans
        small values: that round trip is lossy, and the whole point here is to reproduce the matrix product
        exactly. A tuple (P, Q) of Vector{Float64} over the full grid is returned.
"""
function expandBVector(vec::Vector{Float64}, primitives::Bsplines.Primitives)
    grid = primitives.grid;    nsL = grid.nsL;    nsS = grid.nsS
    P = zeros( grid.NoPoints );    Q = zeros( grid.NoPoints )
    for  k = 1:nsL
        bs = primitives.bsplinesL[k];   add = 1 - bs.lower
        for  r = bs.lower:min(bs.upper, grid.NoPoints)    P[r] = P[r] + vec[k] * bs.bs[r+add]    end
    end
    for  k = 1:nsS
        bs = primitives.bsplinesS[k];   add = 1 - bs.lower
        for  r = bs.lower:min(bs.upper, grid.NoPoints)    Q[r] = Q[r] + vec[nsL+k] * bs.bs[r+add]    end
    end
    return( (P, Q) )
end


"""
`SelfConsistent.screenedProduct(Vk::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
        primitives::Bsplines.Primitives)`  
    ... forms the vector whose i-th entry is  INT B_i(r) f(r) w_r V_k(r) dr,  with f = P on the large block
        and Q on the small one -- i.e. the product of the screened-potential matrix with a coefficient
        vector, WITHOUT ever building that matrix.

        The matrix that InteractionStrength.XL_CoulombKinkAware assembles has entries
        wm[i,k] = INT B_i B_k w_r V_k, a full double loop over roughly 110000 (i,k) pairs which re-expands
        its B-spline arrays inside the inner loop -- and it is then contracted with a single vector.
        Measured 10-Aug-2026 on Be-like C^2+: that assembly was 0.99 s of a 1.32 s gradient, and the
        gradient was ~80% of a rotation iteration.  Since the matrix is symmetric within each block, both
        products a caller needs use the same potential and differ only in the vector.
        A Vector{Float64} of length nsL+nsS is returned.
"""
function screenedProduct(Vk::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
                                primitives::Bsplines.Primitives)
    grid = primitives.grid;    nsL = grid.nsL;    nsS = grid.nsS
    out  = zeros( nsL + nsS )
    for  i = 1:nsL
        bs  = primitives.bsplinesL[i];    add = 1 - bs.lower
        mtp = min( bs.upper, length(Vk), length(P) );    wa = 0.
        for  r = max(2, bs.lower):mtp    wa = wa + bs.bs[r+add] * P[r] * grid.wr[r] * Vk[r]    end
        out[i] = wa
    end
    for  i = 1:nsS
        bs  = primitives.bsplinesS[i];    add = 1 - bs.lower
        mtp = min( bs.upper, length(Vk), length(Q) );    wa = 0.
        for  r = max(2, bs.lower):mtp    wa = wa + bs.bs[r+add] * Q[r] * grid.wr[r] * Vk[r]    end
        out[nsL+i] = wa
    end
    return( out )
end


"""
`SelfConsistent.computeOrbitalGradient(bVectors::Dict{Subshell, Vector{Float64}},
        coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, nucPot::Radial.Potential,
        storage::Dict{String,Array{Float64,2}})`  
    ... the ANALYTIC gradient dE/db_a of the same energy that energyFromBVectors evaluates, for every
        subshell, as a full B-spline coefficient vector. Nothing here is new machinery: the one-electron
        integral is I(a,b) = b_a^T H1 b_b with H1 = Bsplines.setupLocalMatrix, and the Slater integral is
        R^k(abcd) = b_a^T M b_c with M = InteractionStrength.XL_CoulombKinkAware(k, a, orb_b, c, orb_d,
        primitives) -- the matrix-valued overload that already exists for the Fock build. Each slot in which
        a subshell occurs contributes once, using the symmetry R^k(abcd) = R^k(badc) for the second pair.
        Must be checked against gradientByFiniteDifference before being trusted.
        A Dict{Subshell, Vector{Float64}} is returned.
"""
function computeOrbitalGradient(bVectors::Dict{Subshell, Vector{Float64}},
                                       coeffs1p::Array{Coefficient1p,1},
                                       coeffs2p::Array{Coefficient2p,1},
                                       subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                       matrixB::Array{Float64,2})
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives; canonicalize=false)
    end
    grad = Dict{Subshell, Vector{Float64}}()
    for  sh  in  subshells    grad[sh] = zeros(nsL+nsS)    end

    expanded = Dict{Subshell, Tuple{Vector{Float64},Vector{Float64}}}()
    scale    = Dict{Subshell, Float64}()
    for  sh  in  subshells
        (pR, qR) = SelfConsistent.expandBVector(bVectors[sh], primitives)
        np = min(length(pR), length(orbitals[sh].P));    nq = min(length(qR), length(orbitals[sh].Q))
        num = sum( orbitals[sh].P[1:np] .* pR[1:np] ) + sum( orbitals[sh].Q[1:nq] .* qR[1:nq] )
        den = sum( pR[1:np] .* pR[1:np] )             + sum( qR[1:nq] .* qR[1:nq] )
        scale[sh]    = den > 0. ? num/den : 1.0
        expanded[sh] = (orbitals[sh].P, orbitals[sh].Q)
    end

    # The screened potential depends only on the rank and the ORBITAL PAIR (b,d), not on the slot being
    # differentiated, and many angular coefficients share the same triple.  The orbitals are fixed for the
    # whole of one gradient evaluation, so the memo is local to this call and cannot go stale.  This is the
    # same redundancy that dominates the average-level Fock build (3.5x for argon, 9.1x for Th+), and after
    # the average-level field was memoised the rotation became the larger part of EOL: 66% of it for W+.
    vkCache = Dict{Tuple{Int64,Subshell,Subshell}, Vector{Float64}}()

    # one-electron:  d/db_a [ w * b_a^T H1 b_b ]  contributes to BOTH slots
    h1 = Dict{Int64, Array{Float64,2}}()
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        h1[kappa] = Bsplines.setupLocalMatrix(kappa, primitives, nucPot, storage)
    end
    for  cf  in  coeffs1p
        if  cf.a.kappa != cf.b.kappa    continue    end          ## no one-electron element between kappas
        # The SAME chain-rule factor the two-electron part applies, and for the same reason.  The energy's
        # one-electron term is built from the ORBITALS, and an orbital is scale*expand(b) -- so
        # GrantIab(orb_a, orb_b) = scale_a * scale_b * (b_a^T H1 b_b) and the derivative carries that product.
        # Omitting it is harmless whenever a and b share a sign, which is why it went unnoticed: a DIAGONAL
        # term has scale^2 = +1 always, and a kappa whose subshells happen to be canonicalized alike gives +1
        # too.  It bites only on an OFF-DIAGONAL element between two subshells of OPPOSITE sign, and then it
        # flips the sign of that contribution outright.  Measured on a three-layer Be RAS, validating the
        # gradient along guaranteed-tangent directions: every kappa came out exact to 1.0000 except kappa = -2,
        # the single kappa holding an off-diagonal pair -- 2p_3/2 with scale +1 and 3p_3/2 with scale -1 --
        # where the ratio of finite difference to prediction ran -0.10, -0.56, -0.024.
        w  = cf.T * scale[cf.a] * scale[cf.b]
        hh = h1[cf.a.kappa]
        grad[cf.a] = grad[cf.a] + w * (hh * bVectors[cf.b])
        grad[cf.b] = grad[cf.b] + w * (transpose(hh) * bVectors[cf.a])
    end

    # two-electron:  R^k(abcd) = b_a^T M(a,orb_b;c,orb_d) b_c  and, by R^k(abcd) = R^k(badc),
    #                          = b_b^T M(b,orb_a;d,orb_c) b_d
    #
    # MATRIX-FREE (10-Aug-2026).  M is never formed.  Its entries are INT B_i B_k w_r V_L, so
    # (M v)_i = INT B_i(r) f(r) w_r V_L(r) with f the expansion of v -- and M is symmetric within each
    # block, so the two products a coefficient needs share one potential and differ only in the vector.
    # Each subshell's expansion is built ONCE per gradient rather than per coefficient.
    # The expansions must carry the SAME sign convention as the orbitals the energy is built from.
    # Bsplines.generateOrbitalFromVector canonicalizes each orbital so that P > 0 at small r, and does NOT
    # feed that back into bVectors -- so a raw expansion and its own orbital can be oppositely signed. The
    # screened potential Vk below comes from the ORBITALS (b,d) while the contracted vector came from the
    # RAW b-vector (a,c), and an off-diagonal CSF pair carries an ODD power of a correlating orbital's sign,
    # so the mismatch survives instead of cancelling. Measured on the Be RAS step-2 case, where 2p_1/2 and
    # 2p_3/2 are canonicalized against their raw vectors in every iteration: the analytic gradient agreed
    # with a finite difference to five digits while only s-orbitals were involved, then went 2.5x, 19x and
    # finally SIGN-WRONG as the 2p weight grew -- which is what stalled the line search. Same defect, and
    # the same remedy, as the cVector sign-matching in computeTwoElectronV.
    for  cf  in  coeffs2p
        for  (sA, sB, sC, sD)  in  [ (cf.a, cf.b, cf.c, cf.d), (cf.b, cf.a, cf.d, cf.c) ]
            # the same triangular-delta and parity guard XL_CoulombKinkAware applies before doing any work
            lA = Basics.subshell_l(sA);   jA = Basics.subshell_2j(sA)
            lB = Basics.subshell_l(sB);   jB = Basics.subshell_2j(sB)
            lC = Basics.subshell_l(sC);   jC = Basics.subshell_2j(sC)
            lD = Basics.subshell_l(sD);   jD = Basics.subshell_2j(sD)
            if  AngularMomentum.triangularDelta(jA+1, jC+1, cf.nu+cf.nu+1) *
                AngularMomentum.triangularDelta(jB+1, jD+1, cf.nu+cf.nu+1) == 0   ||
                rem(lA+lC+cf.nu, 2) == 1   ||   rem(lB+lD+cf.nu, 2) == 1     continue
            end
            xc = AngularMomentum.CL_reduced_me(sA, cf.nu, sC) * AngularMomentum.CL_reduced_me(sB, cf.nu, sD)
            if  rem(cf.nu, 2) == 1    xc = - xc    end
            Vk = get!(vkCache, (cf.nu, sB, sD)) do
                     RadialIntegrals.buildScreenedPotential(cf.nu, orbitals[sB], orbitals[sD], primitives.grid;
                                                                  mtpOut=primitives.grid.NoPoints)
                 end
            (pC, qC) = expanded[sC];     (pA, qA) = expanded[sA]
            grad[sA] = grad[sA] + (cf.V * xc * scale[sA]) * SelfConsistent.screenedProduct(Vk, pC, qC, primitives)
            grad[sC] = grad[sC] + (cf.V * xc * scale[sC]) * SelfConsistent.screenedProduct(Vk, pA, qA, primitives)
        end
    end

    # THE SCALE-INVARIANCE PROJECTION.  generateOrbitalFromVector NORMALIZES, so the functional depends on b
    # only through b/sqrt(b^T B b) and is invariant under b -> lambda b.  Differentiating that identity at
    # lambda = 1 gives <grad, b> = 0 for the true gradient, and the correction is a projection along B b:
    #     grad  <-  grad - (<grad,b> / b^T B b) * (B b).
    # Along B b, not along b -- the metric matters, and getting it wrong was one of the wrong turns of
    # 31-Aug-2026.  The virtual directions are already B-orthogonal to the occupied orbitals, so on THEM this
    # correction is exactly zero and neither the direction nor gNorm change.  It bites on `dir`, which is a
    # CONJUGATE-GRADIENT combination carrying dirPrev from the previous iteration's basis, built against the
    # previous orbitals: measured b^T B dir / (|b|_B |dir|_B) = +3.9e-02 and -1.2e-02 on the Be RAS step-2
    # case, and that leak is exactly the 3 % by which <grad,dir> exceeded a converged central difference
    # after item 121 was fixed (ratio 0.970348, stable over four decades of step size).
    for  sh  in  subshells
        bb = transpose(bVectors[sh]) * matrixB * bVectors[sh]
        if  abs(bb) > 1.0e-30
            grad[sh] = grad[sh] - ( sum( grad[sh] .* bVectors[sh] ) / bb ) * (matrixB * bVectors[sh])
        end
    end
    return( grad )
end


"""
`SelfConsistent.gradientByFiniteDifference(bVectors::Dict{Subshell, Vector{Float64}},
        virtuals::Dict{Subshell, Array{Vector{Float64},1}}, coeffs1p::Array{Coefficient1p,1},
        coeffs2p::Array{Coefficient2p,1}, subshells::Array{Subshell,1},
        primitives::Bsplines.Primitives, grid::Radial.Grid, nucPot::Radial.Potential; hStep::Float64=1.0e-4)`  
    ... the gradient of the EOL energy with respect to the allowed orbital rotations, by central finite
        differences. Slow but free of any derivation, so it is the reference against which an analytic
        gradient must later be checked, and on its own it answers the question "is this converged solution
        actually stationary?". Each direction is S-orthogonal to the orbital itself, so normalization
        contributes only at second order and is omitted. A Dict{Subshell, Vector{Float64}} is returned.
"""
function gradientByFiniteDifference(bVectors::Dict{Subshell, Vector{Float64}},
                                           virtuals::Dict{Subshell, Array{Vector{Float64},1}},
                                           coeffs1p::Array{Coefficient1p,1},
                                           coeffs2p::Array{Coefficient2p,1},
                                           subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                           grid::Radial.Grid, nucPot::Radial.Potential; hStep::Float64=1.0e-4)
    grad = Dict{Subshell, Vector{Float64}}()
    for  sh  in  subshells
        dirs = virtuals[sh];    g = zeros( length(dirs) )
        for  (iv, phi)  in  enumerate(dirs)
            bPlus  = copy(bVectors);    bPlus[sh]  = bVectors[sh] + hStep * phi
            bMinus = copy(bVectors);    bMinus[sh] = bVectors[sh] - hStep * phi
            ePlus  = SelfConsistent.energyFromBVectors(bPlus,  coeffs1p, coeffs2p, subshells, primitives, grid, nucPot)
            eMinus = SelfConsistent.energyFromBVectors(bMinus, coeffs1p, coeffs2p, subshells, primitives, grid, nucPot)
            g[iv]  = (ePlus - eMinus) / (2 * hStep)
        end
        grad[sh] = g
    end
    return( grad )
end


"""
`SelfConsistent.projectOntoPositiveBranch(bVectors::Dict{Subshell, Vector{Float64}},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, nucPot::Radial.Potential,
        matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}})`  
    ... projects every orbital onto the POSITIVE-ENERGY branch of its own kappa block, then S-orthonormalizes
        the orbitals of each kappa among themselves, entirely in B-spline coefficient space. Returns the
        projected (bVectors, maximum negative-branch weight removed).

        This is what keeps the rotation from collapsing into the Dirac sea. The Dirac operator is unbounded
        below, so the energy has no minimum once an orbital acquires negative-continuum character. The
        rotation DIRECTIONS are all positive-branch by construction, but the round trip
        b -> Bsplines.generateOrbitalFromVector -> b truncates at mtp and cleans small values, and that
        residue carries negative-branch content. Measured for Be with 48 directions: the 2s negative-branch
        weight jumps from 5e-12 to 1.2e-04 in one step and to 3.3e-03 in the next, at which point <h1> for 2s
        has gone from -1.6 to -43.5 and the energy to -69.9 Ha. Projecting after every step removes the
        residue before it can be amplified. The orthonormalization is done here, in coefficient space, rather
        than via orthonormalizeSameKappa, precisely to avoid that lossy round trip.
"""
function projectOntoPositiveBranch(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
                                          primitives::Bsplines.Primitives, nucPot::Radial.Potential,
                                          matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}};
                                          spectrum=nothing, frozen::Array{Subshell,1}=Subshell[])
    out    = Dict{Subshell, Vector{Float64}}()
    worst  = 0.
    posSpec = isnothing(spectrum) ?
              SelfConsistent.positiveBranchSpectrum(subshells, primitives, nucPot, matrixB, storage) : spectrum
    posSet  = Dict{Int64, Array{Vector{Float64},1}}()
    for  kappa  in  unique( [sh.kappa for sh in subshells] )    posSet[kappa] = posSpec[kappa][1]    end
    # (1) project each orbital on the positive branch of its kappa.  A FROZEN orbital is passed through
    # untouched: it arrives from a converged computation, so it is already on the positive branch and already
    # normalized, and projecting it again can only move it.
    frozenSet = Set(frozen)
    for  sh  in  subshells
        if  sh in frozenSet    out[sh] = copy(bVectors[sh]);    continue    end
        b = bVectors[sh];    v = zeros( length(b) )
        for  phi  in  posSet[sh.kappa]    v = v + (transpose(phi) * matrixB * b) * phi    end
        nrm2Full = abs( transpose(b) * matrixB * b );    nrm2Pos = abs( transpose(v) * matrixB * v )
        worst    = max( worst, 1.0 - nrm2Pos/max(nrm2Full,1.0e-30) )
        out[sh]  = v / sqrt( max(nrm2Pos, 1.0e-30) )
    end
    # (2) S-orthonormalize within each kappa, in coefficient space, so the positive span is preserved --
    # but ONLY where it is actually needed.  Gram-Schmidt is sequential and asymmetric: it leaves the first
    # orbital of a kappa untouched and pushes the whole correction onto the later ones, so applying it to
    # an already-orthonormal set rotates the orbitals for nothing.  Measured on Li, doing it unconditionally
    # cost 8.5e-07 Ha while removing a negative-branch weight of only 3.4e-14 -- more than the entire
    # discrepancy that sent us looking.  Skip it when the block is orthonormal to tolerance.
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        shkAll = [ sh for sh in subshells if sh.kappa == kappa ]
        frz    = [ sh for sh in shkAll if      sh in frozenSet ]
        shk    = [ sh for sh in shkAll if  !( sh in frozenSet ) ]
        # THE FROZEN MEMBERS OF THIS BLOCK ARE HELD FIXED, and the actives are made orthogonal to THEM before
        # anything else.  Loewdin below treats every orbital of a block alike -- that is its virtue, and here
        # its problem: applied to the whole block it rotates the frozen orbitals too, and on a five-step RAS
        # expansion at Z = 92 that moved the frozen 1s by 6.2 % at the fifth step, where the kappa = -1 block
        # first holds six s orbitals and the deviation finally exceeds the 1e-9 gate.  The energy of a state
        # still 98 % reference then rose by 568 Ha.  So: project the actives onto the orthogonal complement of
        # the frozen set first, then Loewdin among the ACTIVES only.
        for  sh  in  shk
            for  f  in  frz
                out[sh] = out[sh] - ( transpose(out[f]) * matrixB * out[sh] ) * out[f]
            end
            nrm = sqrt( abs( transpose(out[sh]) * matrixB * out[sh] ) )
            if  nrm > 1.0e-12    out[sh] = out[sh] / nrm    end
        end
        isempty(shk)  &&  continue
        dev = 0.
        for  (i, sha) in enumerate(shk),  (j, shb) in enumerate(shk)
            ov  = transpose(out[sha]) * matrixB * out[shb]
            dev = max( dev, abs( ov - (i == j ? 1.0 : 0.0) ) )
        end
        if  haskey(ENV, "JAC_EOL_LOWDIN")
            @printf(">> [EOL-LOWDIN] kappa %3d  dev %.6e  %s\n", kappa, dev,
                    dev < 1.0e-9 ? "-> SKIPPED (below the 1e-9 gate)" : "-> Loewdin FIRES");   flush(stdout)
        end
        if  dev < 1.0e-9    continue    end
        # SYMMETRIC (Loewdin) orthogonalisation, S^(-1/2), in place of Gram-Schmidt.  Gram-Schmidt is
        # sequential and asymmetric: it leaves the FIRST orbital of a kappa untouched and pushes the whole
        # correction onto the later ones, so perturbing an earlier orbital of a kappa silently moves the later
        # ones as well.  Loewdin treats every orbital of the block alike and is the orthonormal set CLOSEST to
        # the input in a least-squares sense, so it introduces no ordering of its own.
        nk  = length(shk)
        ovl = zeros(nk, nk)
        for  (i, sha) in enumerate(shk),  (j, shb) in enumerate(shk)
            ovl[i,j] = transpose(out[sha]) * matrixB * out[shb]
        end
        ovl = 0.5 * (ovl + transpose(ovl))          ## exact symmetry before the eigendecomposition
        wa  = LinearAlgebra.eigen(ovl)
        if  haskey(ENV, "JAC_EOL_LOWDIN")
            @printf(">> [EOL-LOWDIN] kappa %3d  overlap eigenvalues: %s   (cut at 1e-12, %d dropped)\n",
                    kappa, join([@sprintf("%.4e", x) for x in wa.values], " "),
                    count(x -> x < 1.0e-12, wa.values));    flush(stdout)
        end
        sinv = zeros(nk, nk)
        for  k = 1:nk
            if  wa.values[k] < 1.0e-12    continue    end     ## a linearly dependent block keeps its input
            sinv = sinv + (wa.vectors[:,k] * transpose(wa.vectors[:,k])) / sqrt(wa.values[k])
        end
        newv = [ zeros(length(out[shk[1]]))  for i = 1:nk ]
        for  i = 1:nk,  j = 1:nk    newv[i] = newv[i] + sinv[j,i] * out[shk[j]]    end
        for  (i, sh) in enumerate(shk)
            nrm = sqrt( abs(transpose(newv[i]) * matrixB * newv[i]) )
            if  nrm > 1.0e-10    out[sh] = newv[i] / nrm    end
        end
    end
    return( out, worst )
end


"""
`SelfConsistent.solveOptimizedLevelFieldByRotation(basis::Basis, nuclearModel::Nuclear.Model,
        primitives::Bsplines.Primitives, settings::AsfSettings; printout::Bool=true)`  
    ... EOL by direct minimization of the energy over orbital ROTATIONS, instead of by solving a per-subshell
        generalized eigenvalue problem. It sits beside the older solveOptimizedLevelField, which is left
        untouched; this one is what SelfConsistent.performSCF dispatches to for Basics.EOLField, and hence
        what every layer of a RasExpansion runs.

        Why: the eigenvalue formulation is ill-posed exactly when a correlating CSF's weight becomes small.
        Measured (09-Aug-2026) for Be 1s^2 2s^2 + 1s^2 2p^2, the present scheme converges to a DEGENERATE
        stationary point -- the 2p^2 weight collapses to -4e-5, so the 2p orbital no longer enters the energy
        at all and its gradient vanishes for a trivial reason (|g| = 9e-7). Its energy is 0.0195 Ha above the
        true minimum. Minimizing over rotations cannot fall into that trap: the orbital is driven by the
        energy gradient rather than selected from a spectrum, and it is free to change character -- which a
        correlation orbital must be, and which maximum-overlap selection (tested, refuted) forbids.

        Each outer step re-diagonalizes the CI, rebuilds the angular coefficients from the target level(s),
        and takes one backtracking-line-search step along the negative gradient projected on the
        active-virtual rotations. Occupied-occupied rotations are excluded as redundant.
        A multiplet::Multiplet of the target block(s) is returned.

        settings.frozenSubshells IS HONOURED: the orbitals it names are held exactly fixed and take no part in
        the gradient, the search direction, the curvature history or the line search, while still entering the
        energy and the CI matrix and still being orthogonalized against. That distinction is what a RAS layer
        means by freezing an orbital, and Basics.generate sets the field for every step of a RasExpansion.
        Freezing everything is allowed and is reported as such: the returned multiplet is then the CI result on
        the orbitals as given.

        WHAT "CONVERGED" MEANS HERE, measured 10-Sep-2026 and worth reading before quoting a number to better than
        a milli-Hartree. THE ROUTE FINDS A STATIONARY POINT, NOT NECESSARILY THE MINIMUM, and which one it finds
        depends on where it started. Three starting bases on three systems, one route, one tolerance:

            case      average-level start   Fock-point start    hydrogenic start    spread(AL, Fock)
            Be Z=4      -14.619514867        -14.619561354       -14.580896847        0.047 mHa
            Be Z=26    -812.786662010       -812.787235260      -812.774914173        0.573 mHa
            C  Z=6      -37.722777567        -37.722694227       -37.681361605        0.083 mHa

        ALL NINE REPORTED CONVERGED. Three things follow, none of them a defect to be fixed:

        * A HYDROGENIC START COSTS 12 TO 41 mHa, which is why SelfConsistent.performSCF always starts this route
          from an average-level basis and never from the initial guess. That decision predates this measurement;
          the numbers are what it is worth.
        * BETWEEN THE TWO LEGITIMATE STARTS the answers differ by 0.05 to 0.57 mHa, and NEITHER IS SYSTEMATICALLY
          BETTER -- the Fock start wins on both Be cases, the average-level start on C. That is a flat valley with
          several shallow stationary points, not one start stopping early.
        * TIGHTENING accuracyScf DOES NOTHING. At 1e-8, 1e-10 and 1e-12 the Be Z=26 answers are BIT-IDENTICAL in
          identical iteration counts, because the exit is not accuracyScf at all: it is the energy's own
          double-precision resolution, 32 eps |E| = 3.64e-12 Ha there, four orders of magnitude below the
          tolerance. Each run is converged as far as arithmetic permits and they still land 0.573 mHa apart.

        So CONVERGED means "stationary to machine resolution FROM THIS START", not "the minimum". The variational
        principle still holds -- a lower energy is the better calculation -- but nothing here identifies the lowest
        of these as global. Where a milli-Hartree matters, run from more than one start and keep the lowest; it
        costs one extra run and is the only honest way to bound the answer.
"""
function solveOptimizedLevelFieldByRotation(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                         settings::AsfSettings; printout::Bool=true, nVirtual::Int64=16,
                                         method::Symbol=:lbfgs, cgRestart::Int64=10, lbfgsMemory::Int64=12)
    # Measurement hook only; the default is unchanged.  The code's own note says conjugacy tolerates the
    # per-iteration rebuild of virtualDirections because it carries ONE previous direction and restarts every
    # ten, while L-BFGS accumulates several pairs and does not -- and its y-pairs are differences of gVec,
    # the same metric-mixed object that item 121 had to stop using for the directional derivative.  This lets
    # that be tested without touching the default.
    # THE ROUTE SUPPLIES THE STEP SPACE, THE STEPPING AND THE ITERATION BUDGET.  A settings carrying an
    # AutomaticRoute leaves all three at the values below, so naming only the field remains a complete
    # specification and behaves exactly as it did before the route existed.  The stepping names are Rule 21's,
    # which describe what the method does; :curvature is the short-memory secant method the code calls :lbfgs.
    if  settings.scfRoute isa Basics.RotationRoute
        nVirtual = settings.scfRoute.nVirtual
        method   = settings.scfRoute.stepping == :curvature ? :lbfgs : settings.scfRoute.stepping
    end
    if  haskey(ENV, "JAC_EOL_METHOD")    method = Symbol(ENV["JAC_EOL_METHOD"])    end
    # Measurement hook only; the default is unchanged.  The search direction is assembled inside the span that
    # virtualDirections returns, and that span holds at most nVirtual vectors out of a b-space of nsL + nsS --
    # 16 of 183 for a C-like uranium RAS layer, i.e. 9 %.  gNorm is built from the components along the SAME
    # span, so the reported |grad| is the REACHABLE gradient and not the gradient.  Whether the descent slows
    # because the surface is flat or because the step space is truncated cannot be told apart without varying
    # this, and it has never been varied: the two callers in module-SelfConsistent.jl both take the default.
    if  haskey(ENV, "JAC_EOL_NVIRTUAL")
        nVirtual = something(tryparse(Int, ENV["JAC_EOL_NVIRTUAL"]), nVirtual)
    end
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    nucPot  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)

    bVectors = Dict{Subshell, Vector{Float64}}()
    for  sh  in  basis.subshells
        bVectors[sh] = Bsplines.fitVectorToPrimitives(basis.orbitals[sh], primitives, matrixB)
    end

    if  settings.levelSelectionCI.active  &&  !isempty(settings.levelSelectionCI.symmetries)
        relevantSyms = unique( settings.levelSelectionCI.symmetries )
    else
        relevantSyms = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
    end
    blockCaches = Dict()
    for  sym  in  relevantSyms    blockCaches[sym] = SelfConsistent.cacheCsfPairCoefficientsEOL(sym, basis)   end

    # Built ONCE: it depends only on the nuclear potential and the B-spline basis, never on the orbitals.
    posSpectrum = SelfConsistent.positiveBranchSpectrum(basis.subshells, primitives, nucPot, matrixB, storage)
    # A FROZEN ORBITAL MUST SURVIVE THIS PROJECTION UNCHANGED (item 32, 01-Sep-2026), and the actives must
    # still come out orthogonal to it -- which is why the frozen set is passed IN rather than restored
    # afterwards.  Restoring afterwards was tried first and is wrong: the projection orthogonalizes the whole
    # kappa block jointly, so putting the frozen vectors back leaves the actives orthogonal to vectors that no
    # longer exist, and a non-orthogonal CSF basis can return an arbitrarily low energy.  Measured: it turned
    # the +568 Ha of the five-step case into -5.28 Ha, a "gain" five thousand times the layer before it.
    (bVectors, _) = SelfConsistent.projectOntoPositiveBranch(bVectors, basis.subshells, primitives,
                                                                     nucPot, matrixB, storage; spectrum=posSpectrum,
                                                                     frozen=settings.frozenSubshells)
    # FROZEN SUBSHELLS.  settings.frozenSubshells names orbitals that must NOT be varied, which is what a RAS
    # layer means by "frozen" -- Basics.generate sets it for every step of a RasExpansion, so this is the
    # ordinary case for this driver and not an exotic one.  The optimizer below therefore runs over
    # activeSubshells alone, while the ENERGY, the CI matrix, the virtual directions and the positive-branch
    # projection continue to see EVERY subshell: a frozen orbital still contributes to the energy and still
    # has to be orthogonalized against.  Nothing has to be re-imposed afterwards, because virtualDirections
    # builds each subshell's directions orthogonal to all occupied orbitals of the same kappa, frozen ones
    # included -- so a step taken by an active orbital leaves the kappa block orthonormal by construction.
    # The pinned vectors are the ones RESTORED just above, i.e. exactly as they were handed in, and they are
    # written back after every later projection too: projectOntoPositiveBranch renormalizes, and its
    # Gram-Schmidt rotates whatever it is given on the occasions it fires.  Until 01-Sep-2026 they were taken
    # AFTER the initial projection instead, on the reasoning that they should "sit on the positive branch like
    # everything else" -- which is true of an orbital built here, but not of one inherited from a converged
    # computation that already satisfies both conditions.  See the measurement above.
    frozenSubshells = [ sh  for sh in basis.subshells  if    sh in settings.frozenSubshells ]
    activeSubshells = [ sh  for sh in basis.subshells  if  !(sh in settings.frozenSubshells) ]
    # The descent test and the radial memo both work on the SPLIT functional; see
    # SelfConsistent.computeFunctionalSplit for why.  frozenRk survives the whole run, because a Slater
    # integral over four frozen orbitals cannot change while those orbitals do not.
    frozenSet   = Set(frozenSubshells)
    isFrozenSub = sh -> sh in frozenSet
    frozenRk    = Dict{NTuple{5,Any}, Float64}()
    pinnedB         = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in frozenSubshells )
    restoreFrozen!  = function(bs)
        for  sh  in  frozenSubshells    bs[sh] = copy(pinnedB[sh])    end
        return( bs )
    end
    # A COST WARNING BEFORE THE EXPENSIVE PART, so that a run which cannot finish says so in the first second
    # instead of taking the machine down without writing a line to the log.  Three runs of a 6 163-CSF layer were
    # killed by the OOM reaper at 12.9, 14.1 and 14.7 GB after SEVEN MINUTES OF COMPLETE SILENCE, and all three
    # looked like a dropped network connection rather than a memory failure; that is what this exists to prevent.
    #
    # THE LAW IS MEASURED, NOT GUESSED -- AND IT WAS RE-MEASURED 04-Sep-2026, WHICH CORRECTED IT TWICE OVER.
    # A plain Coulomb CI of the same ion fits 1.58 GB + 0.181 MB/CSF; the EOL field costs about ten times that
    # per CSF because it carries per-level b-vectors forward across iterations.
    #   THE OLD LAW, 1.58 + 1.8e-3 nCsf, PREDICTED 12.7 GB FOR THE OS^16+ CASE THAT DIED AND LOOKED RIGHT
    # AGAINST "12.9-14.7 GB observed" -- but those three numbers were OOM KILLS on a 15 GB machine, i.e. the
    # ceiling and not the requirement.  Run under a 25 GB cap on a machine that has it, the same case needs
    # 18.29 GB and COMPLETES in 1270 s.  It was never diverging.
    #   AND THE COST DEPENDS ON THE NUMBER OF TARGET LEVELS, which the old law did not carry.  Two points,
    # same basis (6163 CSFs), RSS sampled from outside the process once a second:
    #        1 target level    11.13 GB    956 s
    #       13 target levels   18.29 GB   1270 s
    # so 0.60 GB per extra level on an 11.1 GB base.  The suspicion the item was opened on -- that the cost
    # grows with (target levels x iterations) RATHER than with the basis -- is therefore only half right: the
    # levels do cost, and the BASE still dominates.  The two-term law below reproduces both points to 0.02 GB.
    # It is still an ESTIMATE and is labelled as one rather than used to refuse the run, since a refusal on an
    # estimate would be worse than a warning on a good one.
    #
    # THE COST GROWS WITH THE NUMBER OF CSFs AND WITH THE ORBITAL COUNT SEPARATELY, and open d- and f-shells are
    # where both explode at once: the CSF count from the couplings, the orbital work as the square of the subshell
    # count, since the electron-electron repulsion is rebuilt for every PAIR.  Measured on C-like Z = 92 at
    # 4 / 67 / 264 / 658 CSFs over 4 / 9 / 16 / 25 subshells: 21 s, 320 s, 677 s, 2860 s.
    if  printout
        nCsf   = length(basis.csfs);    nSub = length(basis.subshells)
        # THE NUMBER OF TARGET LEVELS ENTERS THE LAW, MEASURED 04-Sep-2026 -- see the note above.  Where the
        # selection is by CONFIGURATION rather than by index the count is not known until the first CI, so 1 is
        # assumed and the estimate is then a LOWER bound; the message says which case it is in.
        nLev   = isempty(settings.levelSelectionCI.indices) ? 1 : length(settings.levelSelectionCI.indices)
        # THE COST IS A LINEAR TERM PLUS A QUADRATIC ONE IN THE SYMMETRY-BLOCK SIZE -- measured 14-Sep-2026,
        # and this is a change of SHAPE, not a re-fit.  The 04-Sep law was 1.58 + 1.0e-3*nCsf*(1.55 +
        # 0.097*(nLev-1)), written per CSF;  fitting that same form to eleven solves over 14-1231 CSFs gave a
        # per-CSF slope 24x smaller.  BOTH ARE RIGHT WHERE THEY WERE FITTED AND WRONG ELSEWHERE, because a
        # single power law cannot describe a sum of two terms:  the linear one dominates below about a thousand
        # CSFs, the quadratic one above, and each fit sampled only its own side.
        #
        # WHY QUADRATIC AT ALL.  Everything built per CSF PAIR within a block -- `PairCoefficientCache`, the CI
        # matrix, the combined-coefficient temporary, the radial caches -- grows with the ordered pairs, i.e. as
        # n^2 for a block of n CSFs.  Item 29 cut the CONSTANT in front of that (the cache is 15x smaller) and
        # could not change the exponent.  The cache is about 38 % of the quadratic term; the other per-pair
        # structures carry the rest, which is why shrinking the cache alone did not move the limit.
        #
        # CALIBRATED ON THREE ANCHORS, Ti III with 3p opened, each run in its own process under /usr/bin/time:
        #
        #        nCsf      sum n_b^2     measured      this law
        #         887       1.74e5        0.15 GB       0.13 GB      (space-dependent part, baseline subtracted)
        #       3 167       2.16e6        0.65 GB       0.64 GB
        #      25 085       1.59e8       24.5   GB     22.1   GB
        #
        # -- within 15 % across a 28x span in CSFs and a 900x span in sum n^2.  At the top anchor the quadratic
        # term is 86 % of the cost;  at the bottom it is 15 %.
        #
        # THE LEVEL COUNT IS DELIBERATELY ABSENT.  The old law charged 0.097 GB per 1000 CSFs per extra target
        # level.  Measured at FIXED space (1231 CSFs, nLev = 1 / 3 / 6): 1.669 / 1.636 / 1.643 GB, i.e. FLAT,
        # where that term predicts +0.6 GB across the range.  Do not reinstate it without a measurement.
        #
        # AND THE BLOCKS ARE REPORTED BECAUSE THEY ARE THE LEVER.  The stores for ALL blocks are held for the
        # whole run, so the cost is the SUM of the squares;  the work, however, is done one block at a time.  On
        # the Ti III balanced row the five blocks hold 6667, 6275, 6092, 4382 and 1669 CSFs, so the largest is
        # 31 % of the sum -- i.e. holding one block at a time would take that row from ~26 GB to ~9 GB.
        blockSizes = Dict{LevelSymmetry,Int64}()
        for  csf  in  basis.csfs
            sym = LevelSymmetry(csf.J, csf.parity);   blockSizes[sym] = get(blockSizes, sym, 0) + 1
        end
        sumN2  = sum( Float64(n)^2  for (_, n) in blockSizes; init=0.0 )
        maxN   = isempty(blockSizes) ? 0 : maximum(values(blockSizes))
        memEol = 1.5 + 1.2e-4 * nCsf + 1.2e-7 * sumN2                    ## GB, measured 14-Sep-2026; see above
        println(">> [EOL-C3] cost estimate: $nCsf CSFs over $nSub subshells, $nLev target level(s)" *
                (isempty(settings.levelSelectionCI.indices) ? " ASSUMED (selection is by configuration, so this is a LOWER bound)" : "") *
                ";  predicted peak " * @sprintf("%.1f GB", memEol) *
                @sprintf(";  %d symmetry blocks, the largest holding %d CSFs and %.0f %% of the pair cost",
                         length(blockSizes), maxN, 100*maxN^2/max(sumN2,1.0)) *
                " -- the stores for ALL blocks are held at once, so that percentage is what one block at a time would cost.")
        # THREADS BUY NOTHING HERE, AND A USER WHO GAVE THEM SHOULD BE TOLD RATHER THAN LEFT TO INFER IT.
        # This solver is single-threaded by construction: it bypasses `Hamiltonian.performCI` -- deliberately,
        # since `diagonalizeBlockEOL` exists to avoid that routine's per-iteration Multiplet-merge overhead --
        # and `performCI` is the only threading JAC has.  A RAS layer is dominated by the EOL rather than by the
        # final CI, so the whole layer runs on one core however many threads were asked for: measured on O III
        # 2p^2 started as `julia -t 6`, %CPU sat at 110, i.e. 1.1 cores of 6.
        #   AND IT IS NOT WORTH REPAIRING, which is why this is a note and not an item.  Threading the
        # per-symmetry block loop was implemented and is CORRECT -- bit-identical to twelve digits at -t 1 and
        # -t 6 -- and WORTHLESS: 62.6 s against 61.5 s, 1.8 %.  The 15-Sep-2026 profile says why: the radial
        # machinery (`buildScreenedPotential`) is the largest identified consumer at 23.3 %, which caps threading
        # at 1.13x on two threads and 1.27x on twelve for the WHOLE run, before any parallel inefficiency.  That
        # is under the standing 30 % bar.  The lever that DOES move this solver is the allocation rate -- 133.5 GB
        # turned over to hold ~2 GB live, with 12 % of wall already in collection -- where removing a single
        # temporary from `combineAngularCoefficientsEOL` bought 1.6x on a 7 062-CSF solve.
        if  Threads.nthreads() > 1  &&  !GBL_EOL_THREAD_NOTE_SHOWN[]
            GBL_EOL_THREAD_NOTE_SHOWN[] = true
            println(">> [EOL-C3] NOTE: Julia was started with $(Threads.nthreads()) threads, and this solver " *
                    "uses ONE.  A RAS layer is dominated by the EOL field, which is single-threaded by " *
                    "construction, so -t N does not speed it up;  measured, %CPU sits at ~110 under -t 6.")
            println(">>   This is measured and deliberate, not an oversight: threading the symmetry-block loop " *
                    "was implemented, verified bit-identical, and gave 1.8 %, and the profile caps any threading " *
                    "of this solver at 1.27x.  Give the threads to a cascade or a process computation instead.")
        end
        if  memEol > 8.0
            println(">> [EOL-C3] *** WARNING: this layer is predicted to need " * @sprintf("%.1f GB", memEol) *
                    ".  A 6 163-CSF layer was killed three times at 12.9-14.7 GB, silently and after seven " *
                    "minutes.  Reduce the layer, the number of target levels, or run where that memory exists.")
            # DELIBERATELY NOT A COLLECTED WARNING.  The warning channel of this function means exactly one
            # thing -- "the EOL field did not converge" -- and the RAS test asserts on precisely that by looking
            # for this function's name in jac-warn.report.  A memory estimate is not a convergence failure, and
            # putting it there would make a large-but-healthy run indistinguishable from a stalled one.
        end
        # AND A TIME ESTIMATE, IN TWO TERMS.  Until 30-Sep-2026 this was a single "4 s/CSF" anchored on C-like U,
        # and a layer's cost is not proportional to its CSF count: there is a FIXED cost per layer, the orbital
        # and pair work, which does not care how many CSFs there are, plus a CI cost which does.  Dividing the
        # total by N_CSF therefore names a quantity that does not exist, and the measured "s/CSF" duly FALLS as
        # the expansion grows -- on Yb+ from 1.045 to 0.184 across four layers.
        #
        #     t  =  A(nSub)  +  b(nSub) * nCsf
        #
        # A = 111 s and b = 0.0426 s/CSF at nSub = 24;  A = 169 s and b = 0.1330 at nSub = 28 (Yb+, 69 electrons,
        # Z = 70, four measured layers).  A rises about as nSub^2.73 and b as nSub^7.39 between those two points.
        # **THOSE TWO EXPONENTS ARE FITTED ON TWO VALUES OF nSub AND ARE NOT THE TRUE SCALING** -- 7.39 is simply
        # what two points give.  They are used because the SHAPE is what matters here and the shape is confirmed
        # out of sample: a fifth Yb+ run at nSub = 33 and 8086 CSFs, both outside the fitted range, came in at
        # 3316 s against 3884 s predicted, 1.17x, where the old law said 42694 s, 12.9x.
        #
        # THE COEFFICIENTS DO NOT TRANSFER BETWEEN ELEMENTS, and that is stated in the printout rather than hidden.
        # Measured on Cl III (15 electrons, Z = 17) at three layers -- 802 CSF / 16 subshells / 114 s, 2816 / 25 /
        # 874 s, 6622 / 36 / 5039 s -- this law gives 0.34x, 0.33x and 1.19x of the measurement, i.e. it UNDER-
        # predicts a light ion by about three where it over-predicts Yb+. Still far better than the old law, which
        # over-predicted those same three layers by 28x, 12.9x and 5.3x.
        #
        # AND THE DIRECTION OF ERROR IS NOW STATED CORRECTLY.  The old message said "within a factor of about two"
        # and "it under-predicts for a larger expansion".  Both were false and in the opposite sense: it OVER-
        # predicted on both systems measured, by 4x to 28x, and the trend is not even monotone in the expansion --
        # worst at the LARGEST layer on Yb+ and at the SMALLEST on Cl III.  A user plans around this number; a
        # 3321-CSF Yb+ layer was billed at 3.7 h, nearly not run for that reason, and took ten minutes.
        # THE EXPONENTS ARE ONLY ANCHORED BETWEEN nSub = 24 AND 28, validated to 33, so they are NOT extrapolated
        # far below that: nSub^7.39 taken down to nSub = 9 is a factor-5000 extrapolation and collapses the CI term
        # to zero -- measured, it billed a 66-CSF / 9-subshell Cl III layer at 8 s against an actual ~25 s, and
        # printed "0.000 s/CSF", which reads as a broken estimator rather than an extrapolated one.  Below the
        # anchored range the law is reported AS an extrapolation and the user is told which way it fails.
        tFixed  = 111.0  * (nSub/24.0)^2.73                   ## the per-layer orbital and pair work
        tPerCsf = 0.0426 * (nSub/24.0)^7.39                   ## the per-CSF CI work
        tEol    = tFixed + tPerCsf * nCsf
        hms(t)  = t < 120. ? @sprintf("~%.0f s", t) : @sprintf("~%.1f h", t/3600)
        println(">> [EOL-C3] time estimate for this layer: " * hms(tEol) * " single-threaded  =  " *
                @sprintf("%.0f s fixed (orbitals, ~nSub^2.7)", tFixed) * " + " *
                @sprintf("%.2e s/CSF x %d CSFs (the CI)", tPerCsf, nCsf) * ".")
        if      nSub < 20
            println(">> [EOL-C3]   *** nSub = $nSub is BELOW the anchored range (24-33), so this is an " *
                    "extrapolation of nSub^7.39 over more than a factor of two and it UNDER-predicts there -- " *
                    "measured 0.3x on a 15-electron ion.  Treat it as a lower bound, not an estimate.")
        elseif  nSub > 40
            println(">> [EOL-C3]   *** nSub = $nSub is ABOVE the anchored range (24-33); the law is untested there.")
        else
            println(">> [EOL-C3]   anchored on Yb+ (69 electrons), exact on its four fitted layers and 1.2x out of " *
                    "sample;  a 15-electron ion measured 0.3x to 1.2x of it, so it is good to a FACTOR OF THREE " *
                    "EITHER WAY across elements.  Its predecessor, one 4 s/CSF anchor, over-predicted by 4x-28x.")
        end
        if  nSub >= 20  &&  nCsf >= 300
            println(">> [EOL-C3] *** NOTE: $nSub subshells and $nCsf CSFs together are the expensive corner -- the " *
                    "repulsion is rebuilt for every orbital PAIR, so the orbital work grows as the square of the " *
                    "subshell count while the CI grows with the CSFs.  For scale, C-like Z = 92 took 2860 s for " *
                    "658 CSFs over 25 subshells, against 21 s for 4 CSFs over 4.")
        end
    end
    if  printout  &&  !isempty(frozenSubshells)
        println(">> [EOL-C3] frozen and NOT optimized: " * join(string.(frozenSubshells), ", ") *
                ";  $(length(activeSubshells)) of $(length(basis.subshells)) subshells are varied.")
    end
    tStep = 1.0;   multiplet = Multiplet("EOL-ByRotation", Level[])
    bestGNorm = Inf;   bestGIter = 0        # the gradient's best, kept for the REPORTED HINT only
    bestE     = Inf;   bestEIter = 0        # the ENERGY's best -- this is what ends the iteration, see below
    # Set by every exit below.  A loop that simply runs out of iterations used to end in silence, which was the
    # fifth of five ways this driver can stop and the only one left unreported.
    stopReason = "";   gNorm = 0.;   iterDone = 0;   eAvailAtExit = NaN;   tStepAtExit = NaN;   eAtExit = NaN
    # THE ACTIVE-PART ENERGY OF EVERY ITERATE, kept so that the layer's DESCENT can be read at any earlier moment.
    # One Float64 per iteration is nothing beside the iterate itself, and it is what turns the doubling test from a
    # second computation into a lookup.  The frozen part is constant within a layer, so a difference of these is a
    # difference of TOTAL energies too.
    eHist = Float64[]
    # THE STEP-HEALTH FLOOR, used by BOTH exits below and hoisted here 05-Sep-2026 so that it can be.
    # A step smaller than this means the line search has collapsed, temporarily or otherwise, and NOTHING
    # measured during such a step says anything about convergence: neither a stagnant gradient nor a flat
    # energy.  Measured (see the stagnation note): a carbon case collapsed to 3.5e-11, recovered a step of
    # 0.125 when allowed to continue, and converged SIXTY-FIVE mHa below where it used to stop.
    stepFloor = 1.0e-6
    e0Prev    = NaN
    collapsedSince = 0        # first iteration of the CURRENT run of collapsed steps; 0 when the step is healthy
    # Direction state, all held in b-space: the virtual directions are rebuilt every iteration, so anything
    # stored in THAT basis would be meaningless one step later.
    dirPrev = Dict{Subshell, Vector{Float64}}();   gPrev = Dict{Subshell, Vector{Float64}}()
    gradPrev = Dict{Subshell, Vector{Float64}}()      # the RAW gradient of the previous iteration; see the
                                                      # curvature pair below for why gVec will not do
    sgPrev  = 0.;    iterSinceRestart = 0
    bPrev   = Dict{Subshell, Vector{Float64}}()
    sHist   = Vector{Dict{Subshell, Vector{Float64}}}();   yHist = Vector{Dict{Subshell, Vector{Float64}}}()
    rhoHist = Float64[]
    # DISCARDING THE CURVATURE WHEN THE CI RE-SOLVE MOVES THE FUNCTIONAL: MEASURED 04-Oct-2026 AND WORSE.
    # Every iteration re-solves the CI eigenvector, so the function the line search minimized at iteration n is
    # not the one it minimizes at n+1, and an L-BFGS pair (s,y) built across that boundary differences the
    # gradients of two DIFFERENT functions.  Measured on 43Ca+ [Ar] 4s with core s -> s singles, that re-solve
    # adds 5.9e-03 Ha of unmodelled descent over 40 iterations, always downward, i.e. 0.2 % to 17 % of the line
    # search's own gain -- so flushing the history whenever it exceeded a fraction of that gain looked obviously
    # right.  IT IS NOT.  At a budget of 150, step 2 of that ladder gives
    #     no flush (as here)        E = -679.517947 Ha        <- best
    #     flush above 10 % of gain  E = -679.517131   (5 flushes,  0.8 mHa worse)
    #     flush above  2 % of gain  E = -679.508011   (19 flushes, 9.9 mHa worse)
    # monotone in the number of flushes: the stale curvature is worth more than the error it carries.  AND THE
    # 2 % RUN IS A TRAP WORTH NAMING -- it returns A = -822.97 MHz against the measured -806.402, the BEST-LOOKING
    # hyperfine constant of the whole study, on the WORST energy of the three.  Judge a solver by its energy.
    # The remedy for the moving functional is to stop alternating at all, i.e. second-order MCSCF; not this.
    for  iter = 1:Basics.maxIterations(settings.scfRoute)
        orbitals = Dict{Subshell, Orbital}()
        for  sh  in  basis.subshells
            orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives; canonicalize=false)
        end
        tempBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
        radial1p  = Dict{Tuple{Subshell,Subshell},Float64}()
        radial2p  = Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}()
        vkCI      = Dict{Tuple{Int64,Subshell,Subshell,Int64},Vector{Float64}}()
        levels    = Level[]
        for  sym  in  relevantSyms
            cache = blockCaches[sym];    idxCsf = cache.idxCsf
            mtx = SelfConsistent.buildCIMatrixEOL(cache, orbitals, grid, nucPot, radial1p, radial2p, vkCI)
            append!( levels, SelfConsistent.diagonalizeBlockEOL(sym, idxCsf, mtx, tempBasis) )
        end
        multiplet    = Basics.sortByEnergy( Multiplet("EOL-ByRotation", levels) )
        targetLevels = SelfConsistent.selectTargetLevelsEOL(multiplet, settings.levelSelectionCI)
        (coeffs1p, coeffs2p) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels)

        # e0 MUST be the energy of the point the line search starts from, measured the way the line search
        # measures its trials -- on the POSITIVE-BRANCH-PROJECTED vectors.  Taking it from the raw bVectors
        # instead made the comparison `eTrial < e0` a comparison between two different functions, and once
        # the two drifted apart no step could ever win: traced on a three-layer Be RAS, the driver reported
        # "no descent found" at iteration 61 while e0 sat 5.42 mHa BELOW the energy of its own starting
        # point, and the trial energy was flat in t across five decades because every t was being measured
        # against that offset.  The direction was never at fault -- dg = -0.047 there, and its tangential
        # part -0.046.  Re-projecting here is a no-op whenever the iterate is already on the manifold, which
        # is what makes it safe: it costs one projection per iteration and removes a whole class of stall.
        (bVectors, _) = SelfConsistent.projectOntoPositiveBranch(bVectors, basis.subshells, primitives,
                                                        nucPot, matrixB, storage; spectrum=posSpectrum)
        restoreFrozen!(bVectors)
        # e0 IS NOW THE ACTIVE PART ALONE.  The frozen part is common to this point and to every trial of the
        # line search below -- the angular coefficients are fixed within this iteration and the frozen orbitals
        # do not move -- so it cancels in `eTrial < e0` exactly, and leaving it out of BOTH sides is what keeps
        # the comparison meaningful at Z = 98, where the total is 3e4 Ha and a step moves its twelfth digit.
        (e0Frozen, e0) = SelfConsistent.energyFromBVectorsSplit(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                         primitives, grid, nucPot, isFrozenSub, frozenRk)
        push!(eHist, e0)
        grad = SelfConsistent.computeOrbitalGradient(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                             primitives, nucPot, storage, matrixB)
        virt = SelfConsistent.virtualDirections(bVectors, basis.subshells, primitives, nucPot,
                                                        matrixB, storage; nVirtual=nVirtual, spectrum=posSpectrum)
        # Project the gradient on the allowed rotations, and PRECONDITION each component by the
        # orbital-energy denominator (eps_v - eps_a) -- the diagonal of the rotation Hessian, i.e. the
        # standard first-order estimate  kappa_av = -g_av / (eps_v - eps_a).  Plain steepest descent
        # converges far too slowly here: it left the Li control 1.8e-6 Ha short after 60 steps.
        # The denominator is floored, since a near-degenerate pair would otherwise produce a huge step.
        gProj = Dict{Subshell, Vector{Float64}}();   step = Dict{Subshell, Vector{Float64}}()
        denom = Dict{Subshell, Vector{Float64}}()
        gNorm = 0.;    sNorm = 0.
        for  sh  in  activeSubshells
            h1k  = Bsplines.setupLocalMatrix(sh.kappa, primitives, nucPot, storage)
            epsA = transpose(bVectors[sh]) * h1k * bVectors[sh]
            gv   = [ transpose(phi) * grad[sh]  for phi in virt[sh] ]
            sv   = zeros( length(gv) );    dv = zeros( length(gv) )
            # MEASURED AND NOT REPAIRED, 02-Sep-2026 -- the preconditioner's denominator is not what its formula
            # assumes, and JAC_EOL_VIRTCHECK below is what shows it.  The division by (epsV - epsA) is the
            # standard orbital-energy denominator, valid when phi is an EIGENSTATE.  It is not: virtualDirections
            # orthogonalizes each direction against the occupied orbitals, so phi^T h1k phi is a Rayleigh quotient
            # of a mixed vector.  On Be-like Z = 92 the sixteen 2s_1/2 directions gave
            #     -5584  -5027  -4457  -3511  -2001  -1993  -935  -360  -164  -131  -110  -54  +438  +635  +734  +3013
            # against an occupied 2s at -1256, where the bare-nucleus s-series is -4861, -1257, -547, -302.  Six lie
            # BELOW the occupied orbital, the lowest beneath even the 1s, and the floor then replaces their large
            # NEGATIVE denominator by the smallest positive one available (+0.05) -- amplifying them and reversing
            # their sign.  Those six carried 100.0 % of the step's weight, while 2p_1/2 and 2p_3/2 had none.  The
            # positive-branch projection then removes them, which is why a planned descent of -35.2 arrives as
            # +3.8e-07 and the search halves 24 times in vain.  Excluding them by  epsV - epsA <= 0  was tried and
            # buys only 1.7e-08 Ha, because the allowed gradient (0.74 of the 1.52 reported) is not reachable
            # through this preconditioner either.  The repair is to build the directions as eigenstates of the
            # mean field so the denominators mean what the formula says; that is a redesign, not a patch.
            ev = zeros( length(gv) )
            for  (iv, phi)  in  enumerate(virt[sh])
                epsV    = transpose(phi) * h1k * phi
                ev[iv]  = epsV
                dv[iv]  = max( epsV - epsA, 0.05 )
                sv[iv]  = - gv[iv] / dv[iv]
            end
            # HOW MUCH OF THE STEP RIDES ON VIRTUALS THAT LIE BELOW THE OCCUPIED ORBITAL?  For such a virtual
            # epsV - epsA < 0, and the floor above replaces that large negative denominator by +0.05 -- the
            # SMALLEST available -- so the component is amplified rather than suppressed, and its sign is
            # reversed.  A negative-energy solution has epsV ~ -2mc^2 ~ -37558 Ha, so this is where a
            # Dirac-sea direction would enter.  Off unless JAC_EOL_VIRTCHECK is set.
            if  haskey(ENV, "JAC_EOL_VIRTCHECK")
                below = [ iv  for iv = 1:length(ev)  if ev[iv] - epsA < 0. ]
                wBelow = isempty(below) ? 0. : sum( sv[below].^2 )
                above  = [ iv  for iv = 1:length(ev)  if ev[iv] - epsA > 0. ]
                gAbove = isempty(above) ? 0. : sqrt(sum( gv[above].^2 ))
                gBelow = isempty(below) ? 0. : sqrt(sum( gv[below].^2 ))
                @printf(">> [EOL-VIRT] iter %d %-9s nVirt = %3d  nBelow = %3d  |g_allowed| = %.6e  |g_below| = %.6e  |s_below|^2/|s|^2 = %.4f\n",
                        iter, string(sh), length(ev), length(below), gAbove, gBelow,
                        sum(sv.^2) > 0. ? wBelow/sum(sv.^2) : 0.)
                flush(stdout)
            end
            gProj[sh] = gv;    step[sh] = sv;    denom[sh] = dv
            gNorm = gNorm + sum( gv.^2 );    sNorm = sNorm + sum( sv.^2 )
        end
        # gNorm SUMS THE GRADIENT OVER THE RETAINED VIRTUAL DIRECTIONS ONLY, so it is the REACHABLE gradient and
        # not the gradient, and two runs with different nVirtual cannot be compared by it -- a larger span
        # mechanically adds terms.  Measured 08-Sep-2026 on C-like uranium layer 2, and recorded here because the
        # symptom invites the opposite conclusion: nVirtual = 48 reports a LARGER |grad| than 16 while reaching a
        # LOWER energy.  The same measurement closed the question of whether |grad| ~ 7.5e-03 is a floor -- IT IS
        # NOT.  That value was simply where a 400-iteration budget stopped; run to 2000 the same case reaches
        # 3.31e-04, with the energy 5.45e-04 Ha (120 cm^-1) deeper.  Enlarging the span does not change the RATE
        # of the descent (16 to 48 identical) but does change where it ENDS (nVirtual 32 is 3.71e-04 Ha lower),
        # and 96 destroys the line search outright at iteration 41.  Do not re-file the plateau as a defect.
        gNorm = sqrt(gNorm);    sNorm = sqrt(sNorm);    iterDone = iter
        # Every subshell frozen is a legitimate request, not an error: the multiplet built at the top of this
        # iteration IS the answer, being the CI result on the given orbitals.  Reported separately because the
        # collapsed-direction exit below would otherwise fire and call it a failure.
        if  isempty(activeSubshells)
            stopReason = "all subshells frozen";   println(">> [EOL-C3] every subshell of this basis is listed in " *
                    "settings.frozenSubshells, so there is nothing to optimize;  the multiplet is the CI result " *
                    "on the orbitals as given.")
            break
        end

        # Assemble the search direction in b-space.  Plain preconditioned steepest descent zigzags here:
        # the energy falls steadily while |grad| merely oscillates, each step undoing part of the previous
        # one.  Polak-Ribiere conjugacy reuses the previous direction to cancel that, at the cost of two
        # dot products and one stored vector per subshell.
        sVec = Dict{Subshell, Vector{Float64}}();   gVec = Dict{Subshell, Vector{Float64}}()
        for  sh  in  activeSubshells
            sv = zeros(nsL+nsS);    gv = zeros(nsL+nsS)
            for  (iv, phi)  in  enumerate(virt[sh])
                sv = sv + step[sh][iv]  * phi
                gv = gv + gProj[sh][iv] * phi
            end
            sVec[sh] = sv;    gVec[sh] = gv
        end
        # THE DEFAULT IS :lbfgs -- Rule 21's CURVATURE stepping -- and it now arrives from the route above.
        # This comment asserted that :conjugate was the default, which it has not been since the five-file
        # split of 27-Aug-2026; corrected 08-Sep-2026 after measuring it, and it was wrong in both halves.
        # On C-like uranium layer 2, :conjugate reaches -14351.681709 against :lbfgs's -14351.757229, i.e.
        # 7.55e-02 Ha WORSE, with |grad| stuck at 3.09e-01 against 7.56e-03.
        # :lbfgs is measurably better where the basis is stable -- Be 1s^2 2s^2 + 1s^2 2p^2 at 12 iterations
        # reaches -14.618710 against conjugacy's -14.616507, i.e. it does change the RATE and not merely the
        # constant -- but it LOSES on the harder Be RAS step-2 case, -14.617374 against -14.619313, and there
        # it discards its own curvature history at iterations 8, 15, 21 and 23.  That is why :conjugate is
        # kept and is worth trying on a case that will not settle.  WHY IT FAILS IS NOT KNOWN.  Two candidates
        # were proposed and BOTH MEASURED SMALL, so neither should be repeated as an explanation:
        #   * the virtual space is NOT churning -- successive frames overlap to |1-<new,old>| = 1e-4..3e-3
        #     with no change in the number of directions, one 0.67 rotation excepted, and that one does not
        #     coincide with any of the four history discards;
        #   * the objective DOES drift, since coeffs1p/coeffs2p are rebuilt from the re-diagonalized mixing
        #     vector every iteration, but only by max|dV| ~ 1e-4..3e-3 against max|V| ~ 1.8.
        # Conjugacy tolerates whatever this is because it carries ONE previous direction and restarts every
        # ten; L-BFGS accumulates five pairs and does not.  Anyone picking this up should measure first.
        dotAll(u, v) = sum( sum(u[sh] .* v[sh])  for sh in activeSubshells )
        # H_0 is the DIAGONAL PRECONDITIONER, not the usual gamma*I: the (eps_v - eps_a) denominators carry
        # real physics and throwing them away for a scalar would be a step backwards.  applyPrecond is what
        # turns gVec into -sVec, so it is exactly the operator already in use.
        applyPrecond = function(q)
            r = Dict{Subshell, Vector{Float64}}()
            for  sh  in  activeSubshells
                rv = zeros(nsL+nsS)
                for  (iv, phi)  in  enumerate(virt[sh])
                    rv = rv + ( (transpose(phi) * q[sh]) / denom[sh][iv] ) * phi
                end
                r[sh] = rv
            end
            return( r )
        end
        # Record the curvature pair (s,y) of the step just taken.  Only pairs with <y,s> > 0 are kept, which
        # is what keeps the implicit inverse Hessian positive definite -- and hence the direction a descent
        # direction -- on a surface that is not convex.
        if  method == :lbfgs  &&  !isempty(bPrev)  &&  !isempty(gradPrev)
            sPair = Dict{Subshell, Vector{Float64}}();   yPair = Dict{Subshell, Vector{Float64}}()
            for  sh  in  activeSubshells
                # y IS A DIFFERENCE OF GRADIENTS, NOT OF gVec.  L-BFGS builds its inverse-Hessian model from
                # (s, y) = (delta x, delta grad f), and gVec = sum_i (phi_i^T grad) phi_i is REBUILT IN A NEW
                # BASIS EVERY ITERATION -- virtualDirections is recomputed against the current orbitals -- so a
                # y taken from it registers curvature whenever the basis merely ROTATES, even where the
                # gradient has not moved.  That spurious curvature is what an accumulated history cannot
                # tolerate and a one-step conjugacy can: measured 31-Aug-2026, L-BFGS died on "no descent" at
                # iteration 179 of Be Scenario B step 3 while :conjugate ran to the limit, though L-BFGS was
                # otherwise the better method (Be A step 2 in 40 iterations against 101, and 1.1 mHa deeper on
                # B before dying).  The raw gradient carries no basis with it.  Same error as item 121, one
                # level down.
                sPair[sh] = bVectors[sh] - bPrev[sh];    yPair[sh] = grad[sh] - gradPrev[sh]
            end
            ys = dotAll(yPair, sPair)
            if  ys > 1.0e-14
                push!(sHist, sPair);    push!(yHist, yPair);    push!(rhoHist, 1.0/ys)
                while  length(sHist) > lbfgsMemory
                    popfirst!(sHist);   popfirst!(yHist);   popfirst!(rhoHist)
                end
            end
        end
        beta = 0.
        dir  = Dict{Subshell, Vector{Float64}}()
        # A DIRECTION MUST NOT TRY TO CHANGE AN ORBITAL'S NORM.  The line search renormalizes every trial
        # vector, so whatever part of the step lies ALONG the orbital is discarded in full before the energy
        # is ever evaluated -- yet `dg = <grad,dir>` counts it, and the Armijo test and the backtracking are
        # then driven by a promise the step cannot keep.  Removing it costs nothing, because it is exactly
        # the part that does nothing.
        #   WHY IT IS NOT ALREADY ZERO.  `virtualDirections` returns a set B-orthogonal to the occupied
        # orbitals, so the STEEPEST-DESCENT direction is clean.  An L-BFGS direction is not: it is a
        # combination of stored (s,y) pairs from earlier iterations, and those were orthogonal to the
        # orbitals as they stood THEN.  The component regrows as the orbitals move, and it regrows fastest
        # for the orbitals that move most -- the correlation virtuals.
        #   MEASURED 04-Oct-2026 on 43Ca+ [Ar] 4s, iteration 20, cos(b, dir) in the B metric:
        #     reference step (no virtuals)   1s..4s, 2p, 3p   5.6e-05 ... 6.1e-04   -- clean
        #     correlation step               1s..4s, 2p, 3p   5.3e-04 ... 7.1e-02
        #                                    5s  0.6155,  6s  0.4023,  7s  0.4459   -- up to 62 % PARALLEL
        # With that much of the step discarded, the realized directional derivative scattered between 0.024
        # and 1.888 of the planned one, against 1.0000 (0.987-1.005) in the reference step; the Armijo ratio
        # measured on the planned step sat at 0.585 while the one measured on the realized displacement sat
        # at 0.824, and it was the realized one that was closer to unity in 114 of 144 accepted steps.  That
        # is the solver's own stated signature for "the model is describing a step the search does not take".
        #   AND IT MUST NOT BE TAKEN FURTHER THAN THE ORBITAL ITSELF -- MEASURED 04-Oct-2026 AND REJECTED.
        # Renormalizing preserves each orbital's NORM but not the mutual ORTHOGONALITY of a kappa block, so the
        # same reasoning seems to say: project the direction against EVERY orbital of its own kappa, not only
        # against the one it moves.  That is WRONG, and the error is physical rather than numerical.  A rotation
        # that mixes 5s into 4s is NOT a redundant degree of freedom while the CSF space is incomplete -- it
        # changes the energy, and it is part of what a correlation layer is for.  Projecting it out over-constrains
        # the optimization: on 43Ca+ [Ar] 4s with core s -> s singles at a budget of 150, step 2 gives
        #     project against the orbital itself only (as here)   E = -679.517947 Ha
        #     project against the whole kappa block               E = -679.498826 Ha   (19 mHa WORSE)
        # and the second returns A = -1009.74 MHz against -806.402 measured, a 25 % error on the worst energy.
        # The component along the orbital ITSELF is the only one the renormalization provably discards, and it is
        # therefore the only one that may be removed here.
        #   AND THE ORTHOGONALITY LOSS THAT MOTIVATED THAT ATTEMPT IS NOT A DEFECT -- SETTLED 04-Oct-2026, so that
        # the 68 % figure is not raised as an alarm again.  A trial vector is b + t*dir renormalized, and
        # virtualDirections has already made dir B-orthogonal to the occupied orbitals of its kappa, so both cross
        # terms in <b1 + t*d1, b2 + t*d2> vanish and the leading surviving term is t^2 <d1,d2>.  The loss is
        # therefore SECOND ORDER IN THE STEP, which is exactly what a first-order method should produce.  Measured
        # on 43Ca+ [Ar] 4s with core s -> s singles, kappa = -1, down one backtracking ladder:
        #     tStep   1.000  0.500  0.250  0.125  0.0625  0.03125  0.015625
        #     dev     0.668  0.397  0.185  0.0669 0.0194  0.00509  0.00129
        # and halving tStep multiplies dev by 0.594, 0.466, 0.361, 0.291, 0.262, 0.253 -> 0.25 = (1/2)^2.  The 0.68
        # belongs to the FIRST trial at tStep = 1.0, which the line search rejects anyway; an ACCEPTED step runs at
        # tStep ~ 0.01, where the deviation is ~1e-4 and Loewdin's correction to it is negligible.
        #   This also explains why that 0.68 came out BIT-IDENTICAL under two different direction rules: at
        # iteration 1 the L-BFGS history is empty, so the direction is the preconditioned steepest-descent step in
        # both, and virtualDirections has already orthogonalized it -- so every variant of this projection is a
        # no-op there and produces the same trial.
        stripNormChange! = function(d::Dict{Subshell, Vector{Float64}})
            for  sh  in  activeSubshells
                bBb = transpose(bVectors[sh]) * matrixB * bVectors[sh]
                if  abs(bBb) > 1.0e-30
                    d[sh] = d[sh] - ( (transpose(bVectors[sh]) * matrixB * d[sh]) / bBb ) * bVectors[sh]
                end
            end
            return( nothing )
        end
        # THE NEWTON DIRECTION -- step (3) of priority item 48, reached with JAC_EOL_METHOD = newton.  It replaces
        # ONLY the choice of direction: the step space, the preconditioner, the positive-branch projection, the
        # frozen-orbital restore, the line search and all five exits stay exactly as the first-order route leaves
        # them, which is what makes the comparison below a comparison of the DIRECTION and of nothing else.
        #
        # WHY A SECOND-ORDER STEP IS THE REMEDY AND NOT A BETTER FIRST-ORDER ONE.  Each iteration re-solves the CI
        # eigenvector, so the function minimized at iteration n is not the one minimized at n+1 and an L-BFGS pair
        # (s,y) differences the gradients of two DIFFERENT functions.  Measured on 43Ca+ [Ar] 4s with core s -> s
        # singles, that re-solve adds 5.9e-03 Ha of unmodelled descent over 40 iterations, always downward.  The
        # cheap repair -- discard the curvature when the relaxation is large -- was implemented and is WORSE,
        # monotonically in the number of flushes.  The REDUCED Hessian contains the relaxation exactly, so the
        # alternation is removed by construction rather than patched around.
        #
        # H_red v = H_oo v  -  2 <x_v| dH/db |c>,   x_v = (H - lambda)^+ P_perp (dH/db . v) c
        #
        # and every piece of it comes from machinery that already exists: H_oo v from a central difference of the
        # EXACT orbital gradient, (dH/db.v)c from a central difference of buildCIMatrixEOL, and the final term from
        # a central difference of the gradient ON THE CI VECTOR -- which is EXACT rather than approximate, because
        # the generalized weight d_rs is quadratic in mc and the gradient linear in the coefficients.  All three
        # routes were verified against each other and against the energy's own mixed second difference before this
        # solver was written (JAC_EOL_HESSOC, step (2)): agreement to 7e-07, and the assembled H_red reproduces the
        # second difference of the RE-SOLVED eigenvalue to ratio 1.000001.
        #
        # ONE H_red v COSTS 4 GRADIENTS AND 2 CI BUILDS.  That is why the inner solver is a TRUNCATED CG and not a
        # formed matrix: the span is nVirtual per active subshell (16 by default), so an explicit Hessian would
        # need one product per column, where CG needs ten to twenty in all.
        #
        # THE SPAN DOES THE REDUNDANCY PROJECTION FOR FREE, which is step (4) of the item and is mostly already
        # answered.  applyPrecond projects onto the span virtualDirections returns, and that span is built
        # B-orthogonal to the occupied orbitals of its own kappa, so every CG vector stays inside it with no extra
        # projection.  What remains is the component along the orbital a direction moves, removed by
        # stripNormChange! below.  Going FURTHER was measured on 04-Oct-2026 and REJECTED -- see the note at
        # stripNormChange!: projecting against the whole kappa block costs 19 mHa on the Ca+ ladder, because a
        # rotation mixing 5s into 4s is NOT redundant while the CSF space is incomplete.
        if      method == :newton  &&  length(targetLevels) != 1
            # The reduced Hessian needs one relaxation term per target level, and two levels of ONE symmetry block
            # couple through it;  neither is built.  Returning a Newton step that silently omits the coupling would
            # be worse than not offering one, so the preconditioned gradient is used and the reason is said once.
            if  iter == 1
                println(">> [EOL-N] JAC_EOL_METHOD = newton is not available for a $(length(targetLevels))-level " *
                        "functional;  the preconditioned first-order step is used instead.  See priority item 48.")
            end
            for  sh  in  activeSubshells    dir[sh] = sVec[sh]    end
        elseif  method == :newton
            lvN     = targetLevels[1]
            symN    = LevelSymmetry(lvN.J, lvN.parity)
            cacheN  = blockCaches[symN];      idxCsfN = cacheN.idxCsf;      nCsfN = length(idxCsfN)
            epsB    = 1.0e-3   ## the b-side difference step, measured as the working value: the re-solved second
                               ## difference degrades to 0.9970 at 1e-5, ordinary round-off for a second difference
            delC    = 1.0e-3   ## the CI-side step;  the difference is EXACT in it, so this has only to avoid round-off
            maxCG   = 20
            cgTol   = 0.1      ## an INEXACT Newton: the caveat of item 48 is that FD quality degrades near
                               ## convergence, and a truncated CG is what tolerates an approximate H*v
            shiftB  = function(v::Dict{Subshell, Vector{Float64}}, t::Float64)
                bb = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
                for  sh  in  activeSubshells    bb[sh] = bVectors[sh] + t * v[sh]    end
                return( bb )
            end
            # The three integral caches MUST be fresh for a perturbed b: they are keyed by bare subshell labels and
            # know nothing of b, so reusing the iteration's caches would return the UNPERTURBED integrals and every
            # difference here would come out zero.
            ciMatAt = function(bv::Dict{Subshell, Vector{Float64}})
                orbs = Dict{Subshell, Orbital}()
                for  sh  in  basis.subshells
                    orbs[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bv[sh], primitives; canonicalize=false)
                end
                return( SelfConsistent.buildCIMatrixEOL(cacheN, orbs, grid, nucPot,
                                Dict{Tuple{Subshell,Subshell},Float64}(),
                                Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}(),
                                Dict{Tuple{Int64,Subshell,Subshell,Int64},Vector{Float64}}()) )
            end
            gradCi  = function(bv::Dict{Subshell, Vector{Float64}}, x::Vector{Float64})
                mcx = zeros( length(basis.csfs) )
                for  (r, idx)  in  enumerate(idxCsfN)    mcx[idx] = x[r]    end
                lvx = Level(lvN.J, lvN.M, lvN.parity, lvN.index, lvN.energy, lvN.relativeOcc,
                            lvN.hasStateRep, lvN.basis, mcx)
                (c1, c2) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, [lvx])
                return( SelfConsistent.computeOrbitalGradient(bv, c1, c2, basis.subshells,
                                                              primitives, nucPot, storage, matrixB) )
            end
            # The target root and its whole block spectrum, re-diagonalized at the CURRENT b so that every quantity
            # below sits at one point.  The spectrum is reused by every CG iteration, so it costs one
            # diagonalization per Newton step and not one per product.  The root is taken by POSITION, since an EOL
            # functional may target an excited root of its block and the position is what survives a perturbation.
            mtxN  = ciMatAt(bVectors)
            eigN  = Basics.diagonalize(MatrixWithLinearAlgebra(), mtxN)
            k0N   = argmin( [ abs(eigN.values[k] - lvN.energy)  for k = 1:length(eigN.values) ] )
            lamN  = eigN.values[k0N]
            evecN = [ eigN.vectors[k] / sqrt(sum(eigN.vectors[k].^2))  for k = 1:nCsfN ]
            cN    = evecN[k0N]
            nHv   = 0
            hessTimes = function(v::Dict{Subshell, Vector{Float64}})
                nHv = nHv + 1
                gp  = gradCi(shiftB(v,  epsB), cN);      gm = gradCi(shiftB(v, -epsB), cN)
                hv  = Dict{Subshell, Vector{Float64}}()
                for  sh  in  activeSubshells    hv[sh] = (gp[sh] - gm[sh]) / (2epsB)    end
                if  nCsfN >= 2
                    # buildCIMatrixEOL fills only the UPPER TRIANGLE and Basics.diagonalize symmetrizes it, so this
                    # product MUST symmetrize too;  with the raw array it drops half of every off-diagonal element
                    # and the first derivative comes out a factor 1.93 wrong, which looks like physics.
                    mp = ciMatAt(shiftB(v, epsB));       mm = ciMatAt(shiftB(v, -epsB))
                    hc = LinearAlgebra.Symmetric((mp - mm) / (2epsB)) * cN
                    rv = hc - sum(cN .* hc) * cN
                    xv = zeros(nCsfN)
                    for  k = 1:nCsfN
                        k == k0N  &&  continue
                        xv = xv + ( sum(evecN[k] .* rv) / (eigN.values[k] - lamN) ) * evecN[k]
                    end
                    # The CI-side difference returns 2 <x|dH/db|c> directly, which is exactly the term wanted.
                    gxp = gradCi(bVectors, cN + delC*xv);    gxm = gradCi(bVectors, cN - delC*xv)
                    for  sh  in  activeSubshells    hv[sh] = hv[sh] - (gxp[sh] - gxm[sh]) / (2delC)    end
                end
                return( hv )
            end
            # PRECONDITIONED, TRUNCATED CG -- AND IT RUNS IN SPAN COORDINATES, NOT IN b-SPACE.  This is the metric
            # trap the note at stripNormChange! below describes for the directional derivative, and it bites a
            # Newton solver twice as hard.  virtualDirections returns a set that is orthonormal in B, NOT in the
            # Euclidean metric, so a CG carried out on b-space vectors with plain dot products silently solves the
            # system in the wrong metric: measured, it left the inner solve hitting all 20 iterations every step and
            # the line search cutting tStep to 0.05, ending 0.28 mHa WORSE than the first-order route.
            #   In span coordinates there is no metric to get wrong.  The unknown is a_i with d = sum_i a_i phi_i;
            # the right-hand side is phi_i^T grad, which is gProj and is already computed above; the operator is the
            # Galerkin projection (H a)_i = phi_i^T (H sum_j a_j phi_j); and the preconditioner is the division by
            # the (eps_v - eps_a) denominators.  Preconditioned steepest descent in these coordinates is EXACTLY
            # -gProj/denom, i.e. the first-order step, so CG starts from what the first-order method would do and
            # can only improve on it if the Hessian model is sound.
            toBSpace = function(a::Dict{Subshell, Vector{Float64}})
                v = Dict{Subshell, Vector{Float64}}()
                for  sh  in  activeSubshells
                    vv = zeros(nsL+nsS)
                    for  (iv, phi)  in  enumerate(virt[sh])    vv = vv + a[sh][iv] * phi    end
                    v[sh] = vv
                end
                return( v )
            end
            applyH = function(a::Dict{Subshell, Vector{Float64}})
                hv = hessTimes( toBSpace(a) )
                r  = Dict{Subshell, Vector{Float64}}()
                for  sh  in  activeSubshells
                    r[sh] = [ transpose(phi) * hv[sh]  for phi in virt[sh] ]
                end
                return( r )
            end
            dotSpan  = function(u, v)
                wa = 0.;    for sh in activeSubshells   wa = wa + sum( u[sh] .* v[sh] )   end;    return( wa )
            end
            precSpan = function(u)
                r = Dict{Subshell, Vector{Float64}}()
                for  sh  in  activeSubshells    r[sh] = u[sh] ./ denom[sh]    end
                return( r )
            end
            aCG = Dict{Subshell, Vector{Float64}}( sh => zeros(length(virt[sh]))  for sh in activeSubshells )
            rCG = Dict{Subshell, Vector{Float64}}( sh => -copy(gProj[sh])         for sh in activeSubshells )
            zCG = precSpan(rCG)
            pCG = Dict{Subshell, Vector{Float64}}( sh => copy(zCG[sh])            for sh in activeSubshells )
            rz  = dotSpan(rCG, zCG);      r0 = sqrt( dotSpan(rCG, rCG) )
            nCG = 0;      negCurv = false;      curv = 0.
            for  it = 1:maxCG
                hp  = applyH(pCG);       pHp = dotSpan(pCG, hp)
                if  it == 1    curv = pHp    end
                # NEGATIVE CURVATURE IS NOT AN ERROR, it is the functional saying the model is not convex here.  CG
                # is truncated at that point and whatever descent it has accumulated is kept;  on the FIRST product
                # there is nothing to keep and the first-order step is used instead.
                if  pHp <= 0.    negCurv = true;    break    end
                al  = rz / pHp
                for  sh  in  activeSubshells
                    aCG[sh] = aCG[sh] + al * pCG[sh];      rCG[sh] = rCG[sh] - al * hp[sh]
                end
                nCG = it
                if  sqrt( dotSpan(rCG, rCG) ) <= cgTol * r0    break    end
                zCG = precSpan(rCG);     rzNew = dotSpan(rCG, zCG)
                for  sh  in  activeSubshells    pCG[sh] = zCG[sh] + (rzNew/rz) * pCG[sh]    end
                rz  = rzNew
            end
            if  nCG == 0
                for  sh  in  activeSubshells    dir[sh] = sVec[sh]    end
            else
                dB = toBSpace(aCG)
                for  sh  in  activeSubshells    dir[sh] = dB[sh]    end
            end
            if  printout
                @printf(">> [EOL-N] iter %d:  CG %d steps, %d H*v products, first curvature %+.4e%s%s\n",
                        iter, nCG, nHv, curv, negCurv ? ", TRUNCATED on negative curvature" : "",
                        nCG == 0 ? ", first-order step used" : "");    flush(stdout)
            end
        elseif  method == :lbfgs  &&  !isempty(sHist)
            # two-loop recursion, giving d = -H grad with H built from the stored pairs around H_0
            q = Dict{Subshell, Vector{Float64}}( sh => copy(grad[sh])  for sh in activeSubshells )
            alphas = zeros( length(sHist) )
            for  i = length(sHist):-1:1
                alphas[i] = rhoHist[i] * dotAll(sHist[i], q)
                for  sh  in  activeSubshells    q[sh] = q[sh] - alphas[i] * yHist[i][sh]    end
            end
            r = applyPrecond(q)
            for  i = 1:length(sHist)
                bb = rhoHist[i] * dotAll(yHist[i], r)
                for  sh  in  activeSubshells    r[sh] = r[sh] + (alphas[i] - bb) * sHist[i][sh]    end
            end
            for  sh  in  activeSubshells    dir[sh] = -r[sh]    end
        elseif  method == :conjugate
            if  !isempty(dirPrev)  &&  iterSinceRestart < cgRestart  &&  abs(sgPrev) > 1.0e-30
                num = 0.
                for  sh  in  activeSubshells    num = num + sum( (gVec[sh] - gPrev[sh]) .* sVec[sh] )    end
                beta = max( 0., num / sgPrev )                   ## Polak-Ribiere+, i.e. restart on beta < 0
            end
            for  sh  in  activeSubshells
                dir[sh] = beta == 0. ? sVec[sh] : sVec[sh] + beta * dirPrev[sh]
            end
        else
            for  sh  in  activeSubshells    dir[sh] = sVec[sh]    end
        end
        # Guard: a conjugate direction must still descend.  If it does not, fall back to steepest descent.
        #
        # THE PAIRING IS WITH THE RAW GRADIENT, NOT WITH gVec, AND THE DIFFERENCE IS A FACTOR OF FIVE TO NINE.
        # virtualDirections returns a B-ORTHONORMAL set (phi^T B phi = 1).  gVec is rebuilt from the gradient's
        # components along that set, gVec = sum_i (phi_i^T grad) phi_i, so pairing it with dir by a plain
        # Euclidean dot gives sum_ij gv_i s_j (phi_i^T phi_j) -- and phi_i^T phi_j is NOT delta_ij, because the
        # phi are orthonormal in B and not in the Euclidean metric.  grad and dir both live in b-space, so
        # their plain pairing IS the directional derivative, with no metric to get wrong.
        #
        # MEASURED 31-Aug-2026 against a central difference of the functional itself (JAC_EOL_FDCHECK), on the
        # Be RAS step-2 case: FD = -6.3452e-03, raw pairing = -6.5391e-03 (3 % high, the residue being the
        # normalization projection that computeOrbitalGradient still omits), old dg = -4.5380e-02 -- SEVEN
        # TIMES the true derivative.  The old value has the right SIGN, which is why the method still descended
        # and why this survived so long; but the line search predicted seven times the decrease it could get,
        # so its first trial always overshot and it halved to 1e-8, which is the plateau seen on Cf^17+ and on
        # Be Scenario B.  gNorm is built from the same components and so never vanished at a stationary point,
        # which is why an energy-based exit test had to exist at all.
        stripNormChange!(dir)
        dg = 0.;   for sh in activeSubshells   dg = dg + sum( grad[sh] .* dir[sh] )   end
        if  dg >= 0.
            for  sh  in  activeSubshells    dir[sh] = sVec[sh]    end
            stripNormChange!(dir)
            beta = 0.
            dg   = 0.;   for sh in activeSubshells   dg = dg + sum( grad[sh] .* dir[sh] )   end
        end
        # ONE-SHOT MEASUREMENT OF A FINITE-DIFFERENCE HESSIAN-VECTOR PRODUCT, off unless JAC_EOL_HESSCHECK is set.
        # The question it answers is whether a second-order (Newton) method needs the orbital-orbital Hessian to be
        # DERIVED at all.  It does not have to be: the gradient here is exact, so H v follows from a central
        # difference of it, (grad(b + eps v) - grad(b - eps v)) / (2 eps), at the cost of two gradients and with no
        # new angular machinery whatever.  Three things are reported, because a Hessian-vector product that is
        # merely plausible is worthless: SYMMETRY, <v1, H v2> against <v2, H v1>, which an exact Hessian must
        # satisfy and a wrong one generally will not; CURVATURE, <v, H v> against the second difference of the
        # ENERGY itself, which ties it to the functional rather than to the gradient routine; and the COST of a
        # gradient against that of an energy, which is what decides whether a Newton-CG inner loop is affordable.
        if  haskey(ENV, "JAC_EOL_HESSCHECK")  &&
                    iter == something(tryparse(Int, ENV["JAC_EOL_HESSCHECK"]), 5)
            hessTimes = function(v::Dict{Subshell, Vector{Float64}}, eps::Float64)
                bp = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
                bm = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
                for  sh  in  activeSubshells
                    bp[sh] = bVectors[sh] + eps * v[sh];    bm[sh] = bVectors[sh] - eps * v[sh]
                end
                gp = SelfConsistent.computeOrbitalGradient(bp, coeffs1p, coeffs2p, basis.subshells,
                                                                   primitives, nucPot, storage, matrixB)
                gm = SelfConsistent.computeOrbitalGradient(bm, coeffs1p, coeffs2p, basis.subshells,
                                                                   primitives, nucPot, storage, matrixB)
                hv = Dict{Subshell, Vector{Float64}}()
                for  sh  in  activeSubshells    hv[sh] = (gp[sh] - gm[sh]) / (2eps)    end
                return( hv )
            end
            dotA = function(x, y)
                wa = 0.;    for sh in activeSubshells   wa = wa + sum( x[sh] .* y[sh] )   end;    return( wa )
            end
            # Two independent probe directions: the search direction and the raw gradient.
            v1 = Dict{Subshell, Vector{Float64}}( sh => copy(dir[sh])   for sh in activeSubshells )
            v2 = Dict{Subshell, Vector{Float64}}( sh => copy(grad[sh])  for sh in activeSubshells )
            for  eps  in  (1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6)
                h1 = hessTimes(v1, eps);    h2 = hessTimes(v2, eps)
                a  = dotA(v2, h1);          b = dotA(v1, h2)
                # <v,Hv> against the energy's own second difference, both on RAW vectors as the gradient is.
                bp = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
                bm = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
                for  sh  in  activeSubshells
                    bp[sh] = bVectors[sh] + eps * v1[sh];    bm[sh] = bVectors[sh] - eps * v1[sh]
                end
                (_, ep) = SelfConsistent.energyFromBVectorsSplit(bp, coeffs1p, coeffs2p, basis.subshells,
                                                     primitives, grid, nucPot, isFrozenSub, frozenRk)
                (_, em) = SelfConsistent.energyFromBVectorsSplit(bm, coeffs1p, coeffs2p, basis.subshells,
                                                     primitives, grid, nucPot, isFrozenSub, frozenRk)
                fd2 = (ep - 2*e0 + em) / eps^2
                @printf(">> [EOL-HESS] eps = %.0e :  sym rel.diff = %.2e ;  <d,Hd> = %+.8e  energy 2nd diff = %+.8e  ratio = %.6f ;  <g,Hg> = %+.6e  curv = %+.6e\n",
                        eps, abs(a-b)/max(abs(a),abs(b),1.0e-30), dotA(v1,h1), fd2,
                        fd2 != 0. ? dotA(v1,h1)/fd2 : NaN, dotA(v2,h2),
                        dotA(v2,v2) > 0. ? dotA(v2,h2)/dotA(v2,v2) : NaN);    flush(stdout)
            end
            tG = time();    for k = 1:5
                SelfConsistent.computeOrbitalGradient(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                              primitives, nucPot, storage, matrixB)
            end;            tG = (time() - tG)/5
            tE = time();    for k = 1:5
                SelfConsistent.energyFromBVectorsSplit(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                     primitives, grid, nucPot, isFrozenSub, frozenRk)
            end;            tE = (time() - tE)/5
            @printf(">> [EOL-HESS] cost: one gradient %.4f s, one energy %.4f s, ratio %.3f;  one H*v = 2 gradients = %.4f s = %.2f energies\n",
                    tG, tE, tG/tE, 2tG, 2tG/tE);    flush(stdout)
        end

        # ONE-SHOT MEASUREMENT OF THE ORBITAL-CI HESSIAN BLOCK AND OF THE REDUCED HESSIAN, off unless
        # JAC_EOL_HESSOC is set.  This is step (2) of the second-order (Newton) MCSCF build order of priority
        # item 48, and it is deliberately a CORRECTNESS GATE placed BEFORE any solver exists: a wrong H_oc would
        # otherwise not announce itself at all, it would look like a convergence problem inside the Newton
        # iteration and be debugged there, which is the expensive place to find it.
        #
        # WHAT IT MEASURES, and why three routes rather than one.  The EOL functional is E(b,c) = c^T H(b) c with c
        # the CI vector, and the alternation that item 48 is about is exactly the OFF-DIAGONAL block H_oc =
        # d2E/db dc: it is the term a first-order method cannot see, because between two iterations the CI re-solve
        # moves c and the gradient history differences two different functions.  H_oc v is obtained here from the
        # EXISTING exact orbital gradient by a central difference ON THE CI VECTOR -- no new angular or radial
        # machinery, exactly as H_oo v needed none (see the HESSCHECK block above).
        #   AND THAT DIFFERENCE IS EXACT, not approximate, which is why no eps-scan is needed for it: the
        # generalized weight d_rs of combineAngularCoefficientsEOL is QUADRATIC in mc and the gradient is LINEAR in
        # the coefficients, so the central difference reproduces the linear term with NO truncation error at all.
        # Only the b-side differences below carry an eps, and they are scanned.
        #
        # THE THREE ROUTES ARE INDEPENDENT CODE PATHS, which is the point:
        #   (a) <v, H_oc w> from the ORBITAL GRADIENT, differenced on c;
        #   (b) the same number as 2 w^T (dH/db . v) c, from buildCIMatrixEOL, differenced on b -- the other
        #       operator order, through the CI matrix rather than through the gradient;
        #   (c) the same number as the MIXED SECOND DIFFERENCE OF THE ENERGY ITSELF, which is what ties it to the
        #       functional rather than to any derivative routine.
        # Agreement of (a) with (b) tests the two implementations against each other; agreement with (c) tests both
        # against the thing being differentiated.  A Hessian block that is merely plausible is worthless.
        #
        # THEN THE ASSEMBLED REDUCED HESSIAN, which is the quantity a Newton step actually needs.  Because c is
        # re-solved at every b, the curvature of the function the solver really walks on is not H_oo but
        #       H_red(v,v) = <v, H_oo v>  -  2 r^T (H - lambda)^+ r ,     r = P_perp (dH/db . v) c
        # the second term being ordinary eigenvalue perturbation theory and NEGATIVE for the lowest root -- which is
        # the curvature-level statement of item 48's measured "always downward" unmodelled descent.  It is checked
        # against the second difference of the RE-SOLVED eigenvalue, i.e. against the true total curvature, which no
        # part of the construction has been told about.
        #
        # EVERYTHING IS RE-DIAGONALIZED AT THE CURRENT bVectors FIRST, and that is not pedantry: `orbitals`,
        # `targetLevels` and the loop's `coeffs1p/coeffs2p` were built BEFORE projectOntoPositiveBranch touched
        # bVectors above.  The projection is a no-op on an iterate already on the manifold, so the two agree in
        # practice -- but a gate that compares a Hessian against an energy must have both at ONE point, or a
        # disagreement cannot be attributed.
        if  haskey(ENV, "JAC_EOL_HESSOC")  &&
                    iter == something(tryparse(Int, ENV["JAC_EOL_HESSOC"]), 5)
            if  length(targetLevels) != 1
                println(">> [EOL-HOC] SKIPPED: this probe is written for ONE target level and the functional has " *
                        "$(length(targetLevels)).  The weighted multi-level case needs one relaxation term per level " *
                        "and the levels of one symmetry block couple; that is a separate measurement.");   flush(stdout)
            else
                lvT    = targetLevels[1]
                symT   = LevelSymmetry(lvT.J, lvT.parity)
                cacheT = blockCaches[symT];     idxCsfT = cacheT.idxCsf;    nCsf = length(idxCsfT)
                dotA   = function(x, y)
                    wa = 0.;   for sh in activeSubshells   wa = wa + sum( x[sh] .* y[sh] )   end;   return( wa )
                end
                # The b-vectors shifted along an orbital direction;  the frozen subshells never move.
                shiftedB = function(v::Dict{Subshell, Vector{Float64}}, t::Float64)
                    bb = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
                    for  sh  in  activeSubshells    bb[sh] = bVectors[sh] + t * v[sh]    end
                    return( bb )
                end
                # The CI matrix of the target block at a given set of b-vectors.  THE THREE INTEGRAL CACHES MUST BE
                # FRESH: they are keyed by bare subshell labels and know nothing of b, so passing the loop's caches
                # would silently return the UNPERTURBED radial integrals and every difference below would be zero.
                ciMatrixAt = function(bv::Dict{Subshell, Vector{Float64}})
                    orbs = Dict{Subshell, Orbital}()
                    for  sh  in  basis.subshells
                        orbs[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bv[sh], primitives; canonicalize=false)
                    end
                    return( SelfConsistent.buildCIMatrixEOL(cacheT, orbs, grid, nucPot,
                                    Dict{Tuple{Subshell,Subshell},Float64}(),
                                    Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}(),
                                    Dict{Tuple{Int64,Subshell,Subshell,Int64},Vector{Float64}}()) )
                end
                # The angular coefficients belonging to an arbitrary block-local CI vector, through the ordinary
                # production path, so the probe measures the code the solver would use and not a copy of it.
                coeffsFromCi = function(x::Vector{Float64})
                    mcx = zeros( length(basis.csfs) )
                    for  (r, idx)  in  enumerate(idxCsfT)    mcx[idx] = x[r]    end
                    lvx = Level(lvT.J, lvT.M, lvT.parity, lvT.index, lvT.energy, lvT.relativeOcc,
                                lvT.hasStateRep, lvT.basis, mcx)
                    return( SelfConsistent.combineAngularCoefficientsEOL(blockCaches, [lvx]) )
                end
                gradFromCi = function(bv::Dict{Subshell, Vector{Float64}}, x::Vector{Float64})
                    (c1, c2) = coeffsFromCi(x)
                    return( SelfConsistent.computeOrbitalGradient(bv, c1, c2, basis.subshells,
                                                                  primitives, nucPot, storage, matrixB) )
                end
                energyFromCi = function(bv::Dict{Subshell, Vector{Float64}}, x::Vector{Float64})
                    (c1, c2) = coeffsFromCi(x)
                    (ef, ea) = SelfConsistent.energyFromBVectorsSplit(bv, c1, c2, basis.subshells,
                                                      primitives, grid, nucPot, isFrozenSub, frozenRk)
                    return( ef + ea )
                end
                # Re-diagonalize HERE, and take the root by its position rather than by "the lowest": an EOL
                # functional may target an excited root of its block, and the position is what stays meaningful
                # when b is perturbed.
                mtx0  = ciMatrixAt(bVectors)
                eig0  = Basics.diagonalize(MatrixWithLinearAlgebra(), mtx0)
                k0    = argmin( [ abs(eig0.values[k] - lvT.energy)  for k = 1:length(eig0.values) ] )
                cVec  = eig0.vectors[k0] / sqrt(sum(eig0.vectors[k0].^2))
                lam   = eig0.values[k0]
                if  nCsf < 2
                    println(">> [EOL-HOC] SKIPPED: the target block holds $nCsf CSF, so there is no CI freedom " *
                            "and H_oc is empty.  Run this on a correlation layer.");   flush(stdout)
                else
                # A CI direction orthogonal to c, deterministic so that a surprising number can be reproduced
                # exactly.  w is normalized;  with |w| = 1 the O(delta^2) norm error of the unnormalized weight
                # d_rs is 1e-06 at delta = 1e-03 and cancels from the antisymmetric differences anyway.
                wVec = [ cos(0.7*r) + 0.3*sin(1.9*r)  for r = 1:nCsf ]
                wVec = wVec - (sum(wVec .* cVec)) * cVec
                wVec = wVec / sqrt(sum(wVec.^2))
                @printf(">> [EOL-HOC] iteration %d:  target %s root %d of %d,  lambda = %+.10f Ha,  |c_target| = %.6f\n",
                        iter, string(symT), k0, nCsf, lam, sqrt(sum(cVec.^2)));   flush(stdout)
                # THE B-METRIC OVERLAP OF THE DIRECTION WITH THE ORBITAL IT MOVES IS PRINTED, and it decides whether
                # the two paths are even comparable.  The CI matrix is built through generateOrbitalFromVector, which
                # NORMALIZES, so H(b) is invariant under b -> lambda b;  the energy and gradient path contracts the
                # coefficients with the b-vector AS GIVEN.  A direction with a component along b would therefore mean
                # two DIFFERENT derivatives and the routes below could not be compared at all.  Measured here it is
                # 1e-16 -- virtualDirections already builds the span B-orthogonal to the occupied orbitals -- so the
                # question is moot and ONE direction suffices.  It is printed rather than assumed because a direction
                # that acquired such a component would make every ratio below wrong for a reason nothing else shows.
                vRaw = Dict{Subshell, Vector{Float64}}( sh => copy(dir[sh])  for sh in activeSubshells )
                nrmV = sqrt( dotA(vRaw, vRaw) );    for sh in activeSubshells   vRaw[sh] = vRaw[sh] / nrmV   end
                for  sh  in  activeSubshells
                    bBv = transpose(bVectors[sh]) * matrixB * vRaw[sh]
                    nb  = sqrt(abs(transpose(bVectors[sh]) * matrixB * bVectors[sh]))
                    nv  = sqrt(abs(transpose(vRaw[sh]) * matrixB * vRaw[sh]))
                    @printf(">> [EOL-HOC]   %-9s b^T B v / (|b|_B |v|_B) = %+.3e   (raw search direction)\n",
                            string(sh), (nb*nv) > 0. ? bBv/(nb*nv) : NaN)
                end
                flush(stdout)
                # BEFORE ANY SECOND DERIVATIVE, THE ZEROTH AND FIRST MUST BE SHOWN TO BE THE SAME FUNCTION on both
                # paths, because every second-derivative route below inherits whichever Hamiltonian its path uses.
                # buildCIMatrixEOL assembles the block from kink-aware XL_CoulombKinkAware integrals; the energy and
                # gradient contract the SAME angular coefficients with the screened-potential route.  If these two
                # differ, a mixed second derivative taken through one and checked against the other must disagree,
                # and no amount of care about eps will reveal which is at fault.
                # IS (values[k0], vectors[k0]) ACTUALLY AN EIGENPAIR?  Asked because it must be, and checked because
                # the whole probe rests on it: every route below pairs that vector with that eigenvalue.
                # buildCIMatrixEOL FILLS ONLY THE UPPER TRIANGLE (`for s = r:n`) and Basics.diagonalize symmetrizes
                # it with LinearAlgebra.Symmetric, so EVERY matrix-vector product here must symmetrize too.  Using the
                # raw array instead drops half of each off-diagonal element: measured on this case it put the first
                # derivative out by a factor 1.93 and left ||H c - lambda c|| = 1.1e-01 on an exact eigenpair.
                hSym  = LinearAlgebra.Symmetric(mtx0)
                resid = sqrt(sum( (hSym*cVec - lam*cVec).^2 ))
                @printf(">> [EOL-HOC] EIGENPAIR :  ||H c - lambda c|| = %.3e ;  c^T H c = %+.10f ;  lambda = %+.10f ;  spectrum = %s\n",
                        resid, sum(cVec .* (hSym*cVec)), lam,
                        string([ @sprintf("%.6f", eig0.values[k])  for k = 1:min(nCsf,5) ]));    flush(stdout)
                ePath = energyFromCi(bVectors, cVec)
                gPath = gradFromCi(bVectors, cVec)
                mP0   = ciMatrixAt(shiftedB(vRaw, 1.0e-4));    mM0 = ciMatrixAt(shiftedB(vRaw, -1.0e-4))
                hco0  = LinearAlgebra.Symmetric((mP0 - mM0) / (2*1.0e-4)) * cVec
                @printf(">> [EOL-HOC] ZEROTH ORDER:  lambda = c^T H_CI c = %+.10f Ha   vs   energy path = %+.10f Ha   difference = %+.3e\n",
                        lam, ePath, lam - ePath)
                @printf(">> [EOL-HOC] FIRST ORDER :  <v,grad> = %+.10e   vs   c^T (dH_CI/db . v) c = %+.10e   ratio = %.6f\n",
                        dotA(vRaw, gPath), sum(cVec .* hco0),
                        sum(cVec .* hco0) != 0. ? dotA(vRaw, gPath)/sum(cVec .* hco0) : NaN);    flush(stdout)
                # AND IF THE FIRST DERIVATIVES DISAGREE, SCAN BOTH FUNCTIONS ALONG THE DIRECTION, because a ratio
                # at one point says only THAT they differ.  Both are evaluated at FIXED c, so this compares the two
                # b-parametrizations alone: f1 is what the line search minimizes, f2 is the expectation value of the
                # matrix whose lowest root the solver reports as the energy.  They agree at t = 0 by construction.
                for  t  in  (0.0, 1.0e-3, -1.0e-3, 1.0e-2, -1.0e-2)
                    bt = shiftedB(vRaw, t)
                    f1 = energyFromCi(bt, cVec)
                    f2 = sum( cVec .* (LinearAlgebra.Symmetric(ciMatrixAt(bt)) * cVec) )
                    @printf(">> [EOL-HOC] SCAN t = %+.0e :  energy path = %+.12f   c^T H_CI c = %+.12f   difference = %+.4e\n",
                            t, f1, f2, f1 - f2);    flush(stdout)
                end
                vDir = vRaw
                # (a) H_oc w from the orbital gradient, differenced on the CI vector.  Exact, as argued above.
                del  = 1.0e-3
                gCiP = gradFromCi(bVectors, cVec + del*wVec);    gCiM = gradFromCi(bVectors, cVec - del*wVec)
                hoc  = Dict{Subshell, Vector{Float64}}()
                for  sh  in  activeSubshells    hoc[sh] = (gCiP[sh] - gCiM[sh]) / (2del)    end
                lhs  = dotA(vDir, hoc)
                # H_oo v v at FIXED CI, for the reduced Hessian below;  this is the quantity step (1) verified.
                gFixP = gradFromCi(shiftedB(vDir,  1.0e-4), cVec)
                gFixM = gradFromCi(shiftedB(vDir, -1.0e-4), cVec)
                hooVV = 0.
                for  sh  in  activeSubshells
                    hooVV = hooVV + sum( vDir[sh] .* (gFixP[sh] - gFixM[sh]) ) / (2*1.0e-4)
                end
                for  eps  in  (1.0e-3, 1.0e-4, 1.0e-5)
                    # (b) the same inner product through the CI matrix, as 2 w^T (dH/db . v) c.
                    mP   = ciMatrixAt(shiftedB(vDir, eps));     mM = ciMatrixAt(shiftedB(vDir, -eps))
                    hco  = LinearAlgebra.Symmetric((mP - mM) / (2eps)) * cVec
                    rhs  = 2 * sum( wVec .* hco )
                    # (c) and as the mixed second difference of the energy, four evaluations.
                    mixed = ( energyFromCi(shiftedB(vDir, eps),  cVec + del*wVec)
                            - energyFromCi(shiftedB(vDir, eps),  cVec - del*wVec)
                            - energyFromCi(shiftedB(vDir, -eps), cVec + del*wVec)
                            + energyFromCi(shiftedB(vDir, -eps), cVec - del*wVec) ) / (4*eps*del)
                    # The CI relaxation of the reduced Hessian, by eigenvalue perturbation theory on the
                    # re-solved root.  The target root is EXCLUDED from the sum, not merely the lowest one.
                    rVec = hco - sum(cVec .* hco) * cVec
                    xVec = zeros(nCsf)
                    for  k = 1:nCsf
                        k == k0  &&  continue
                        ek = eig0.vectors[k] / sqrt(sum(eig0.vectors[k].^2))
                        xVec = xVec + ( sum(ek .* rVec) / (eig0.values[k] - lam) ) * ek
                    end
                    relax = -2 * sum( rVec .* xVec )
                    # The true total curvature: the second difference of the RE-SOLVED eigenvalue.
                    ePls = Basics.diagonalize(MatrixWithLinearAlgebra(), mP).values[k0]
                    eMns = Basics.diagonalize(MatrixWithLinearAlgebra(), mM).values[k0]
                    fdRes = (ePls - 2*lam + eMns) / eps^2
                    @printf(">> [EOL-HOC] eps = %.0e :  <v,H_oc w>  grad-route = %+.8e  CI-route = %+.8e  energy-route = %+.8e  rel.diff = %.2e / %.2e\n",
                            eps, lhs, rhs, mixed, abs(lhs-rhs)/max(abs(lhs),abs(rhs),1.0e-30),
                            abs(lhs-mixed)/max(abs(lhs),abs(mixed),1.0e-30));    flush(stdout)
                    @printf(">> [EOL-HOC]            H_oo = %+.8e   CI relaxation = %+.8e (%.1f %% of H_oo)   H_red = %+.8e   re-solved 2nd diff = %+.8e   ratio = %.6f\n",
                            hooVV, relax, abs(hooVV) > 0. ? 100*abs(relax)/abs(hooVV) : NaN,
                            hooVV + relax, fdRes, fdRes != 0. ? (hooVV + relax)/fdRes : NaN);    flush(stdout)
                end
                end
            end
        end

        # ONE-SHOT FINITE-DIFFERENCE CHECK OF THE GRADIENT, off unless JAC_EOL_FDCHECK is set.
        # Four inferences about this solver's plateau were refuted by measurement on 30/31-Aug-2026, so this
        # measures the thing itself: is <grad,dir> the directional derivative of the functional the line
        # search evaluates?  Both conventions are tested, because the gradient may be that of the raw
        # functional or that of the functional restricted to the normalization manifold, and the line search
        # renormalizes every trial vector.  If a central difference reproduces dg for either convention, the
        # gradient is exact and the tiny steps are the surface's own curvature (conditioning, not a defect);
        # if both differ from dg by a fixed factor, the direction is built on a gradient that is not the
        # functional's, which would explain every plateau seen so far.
        if  haskey(ENV, "JAC_EOL_FDCHECK")  &&
                    iter == something(tryparse(Int, ENV["JAC_EOL_FDCHECK"]), 5)   ## the value selects the iteration
            @printf(">> [EOL-FD] iteration %d:  <grad,dir> = %+.10e\n", iter, dg)
            # IS THE GRADIENT ORTHOGONAL TO ITS OWN b-VECTOR?  generateOrbitalFromVector normalizes, so the
            # functional is SCALE-INVARIANT in b: E(lambda*b) = E(b), and differentiating at lambda = 1 gives
            # <grad, b> = 0 for the true gradient.  A non-zero overlap is a spurious radial component that
            # does nothing to the energy but inflates every directional derivative built from it.  dgProj is
            # the same directional derivative with that component removed, subshell by subshell.
            dgProj = 0.
            for  sh  in  activeSubshells
                bb  = sum( bVectors[sh] .* bVectors[sh] )
                ov  = bb > 0. ? sum( gVec[sh] .* bVectors[sh] ) / bb : 0.
                gp  = gVec[sh] - ov * bVectors[sh]
                cosang = sqrt(sum(gVec[sh].^2)*bb) > 0. ?
                         sum( gVec[sh] .* bVectors[sh] ) / sqrt(sum(gVec[sh].^2)*bb) : 0.
                @printf(">> [EOL-FD]   %-9s cos(grad,b) = %+.6f\n", string(sh), cosang)
                dgProj = dgProj + sum( gp .* dir[sh] )
            end
            @printf(">> [EOL-FD] iteration %d:  <grad_projected,dir> = %+.10e   (ratio to raw %+.6f)\n",
                    iter, dgProj, dg != 0. ? dgProj/dg : NaN)
            # THE CANDIDATE.  virtualDirections returns a B-ORTHONORMAL set (phi^T B phi = 1), but gVec is
            # rebuilt as sum_i (phi_i^T grad) phi_i and then paired with dir by a plain Euclidean dot.  For
            # dir = sum_j s_j phi_j the true directional derivative is sum_i gv_i s_i, whereas that dot gives
            # sum_ij gv_i s_j (phi_i^T phi_j), and phi_i^T phi_j is NOT delta_ij in the Euclidean metric.
            # The raw gradient and the direction both live in b-space, so their plain pairing IS the
            # derivative -- that is what dgRaw measures here.
            dgRaw = 0.
            for  sh  in  activeSubshells    dgRaw = dgRaw + sum( grad[sh] .* dir[sh] )    end
            # IS THE SEARCH DIRECTION B-ORTHOGONAL TO THE ORBITAL IT MOVES?  virtualDirections orthogonalizes
            # every phi against the occupied orbitals of that kappa IN THE B METRIC, so b^T B dir should
            # vanish.  If it does, the scale-invariance correction -- which is a projection along B b, not
            # along b -- contributes nothing to <grad,dir>, and cannot be the 3 % residual against the finite
            # difference.  Printed rather than assumed, because assuming it is how six earlier hypotheses died.
            for  sh  in  activeSubshells
                bBd  = transpose(bVectors[sh]) * matrixB * dir[sh]
                nb   = sqrt(abs(transpose(bVectors[sh]) * matrixB * bVectors[sh]))
                nd   = sqrt(abs(transpose(dir[sh]) * matrixB * dir[sh]))
                @printf(">> [EOL-FD]   %-9s b^T B dir / (|b|_B |dir|_B) = %+.3e\n", string(sh),
                        (nb*nd) > 0. ? bBd/(nb*nd) : NaN)
            end
            @printf(">> [EOL-FD] iteration %d:  <grad_raw,dir> = %+.10e   (ratio to dg %+.6f)\n",
                    iter, dgRaw, dg != 0. ? dgRaw/dg : NaN)
            for  eps  in  (1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7)
                bR = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in basis.subshells )
                bN = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in basis.subshells )
                eR = zeros(2);   eN = zeros(2)
                for  (k, sgn)  in  enumerate((+1.0, -1.0))
                    for  sh  in  activeSubshells
                        v      = bVectors[sh] + sgn*eps * dir[sh]
                        bR[sh] = v
                        bN[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
                    end
                    (_, eR[k]) = SelfConsistent.energyFromBVectorsSplit(bR, coeffs1p, coeffs2p, basis.subshells,
                                                     primitives, grid, nucPot, isFrozenSub, frozenRk)
                    (_, eN[k]) = SelfConsistent.energyFromBVectorsSplit(bN, coeffs1p, coeffs2p, basis.subshells,
                                                     primitives, grid, nucPot, isFrozenSub, frozenRk)
                end
                fdR = (eR[1] - eR[2]) / (2eps);    fdN = (eN[1] - eN[2]) / (2eps)
                @printf(">> [EOL-FD]   eps = %.0e :  raw FD = %+.10e (ratio %+.6f) ;  normalized FD = %+.10e (ratio %+.6f)\n",
                        eps, fdR, dg != 0. ? fdR/dg : NaN, fdN, dg != 0. ? fdN/dg : NaN)
            end
        end
        iterSinceRestart = beta == 0. ? 0 : iterSinceRestart + 1
        sgPrev = 0.;   for sh in activeSubshells   sgPrev = sgPrev + sum( gVec[sh] .* sVec[sh] )   end
        gPrev  = gVec;    gradPrev = grad;    dirPrev = dir
        bPrev  = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in activeSubshells )
        if  sNorm < 1.0e-14
            # Until 18-Aug-2026 this was the one exit of four that said NOTHING, so a run could end here and
            # be read as a completed optimisation.  It means the PRECONDITIONED direction has collapsed, which
            # is not the same as a converged gradient and must not be reported as one.
            stopReason = "direction collapsed";   println(">> [EOL-C3] STOPPED at iteration $iter: the preconditioned direction has collapsed " *
                    "(|s| = $sNorm), with |grad| = $gNorm.  This is NOT convergence.")
            Defaults.warn(AddWarning(), "SelfConsistent.solveOptimizedLevelFieldByRotation(): the EOL field did NOT " *
                          "converge -- the preconditioned direction collapsed at iteration $iter with |grad| = " *
                          @sprintf("%.1e", gNorm) * ".  The energies are NOT self-consistent.")
            break
        end
        if  printout
            println(">> [EOL-C3] iter $iter:  E = $(multiplet.levels[1].energy)   |grad| = $gNorm   step = $tStep")
        end
        # Converged when the GRADIENT is small. settings.accuracyScf is its tolerance, which is the honest
        # test: the old sole criterion, |E - E_prev| < accuracyScf, halts when PROGRESS is slow, which is a
        # statement about the optimizer and not about the solution -- and it is why every EOL value quoted
        # before 16-Aug-2026 was an upper bound.
        #
        # "BOTH TESTS ARE KEPT" IS NO LONGER TRUE, and this sentence stood here after it stopped being true --
        # corrected 04-Sep-2026.  The energy test was REMOVED on 03-Sep in favour of the gradient-STAGNATION
        # exit below (see its own note: a stationary energy stopped every step of Be Scenario A short, at
        # iterations 3, 15, 37 and 80, where the same runs converged at 63, 154, 51 and 40).  THIS IS NOW THE
        # ONLY CONVERGENCE EXIT; every other one below is a give-up and says so.
        #
        # AND IT IS KNOWN TO BE INADEQUATE AT BOTH ENDS OF Z -- priority item 6, measured 03-Sep-2026.  |grad|
        # is not scale-free: its size is set by the excitation energies that set the curvature, ~1000 Ha at
        # Z = 92 against ~1 Ha at Z = 4.  At Z = 92 the run is DONE and cannot say so (energy flat to 4e-15 Ha,
        # a central difference of the functional returning pure round-off, ~3e-09 Ha still available) while
        # |grad| = 0.0023 reads as unconverged against 1e-6; at Z = 4 the exits fire while the energy is still
        # falling by 1.3e-05 Ha per iteration.  The remedy is NOT the naive energy test that was removed -- that
        # was measured and is wrong -- but a GUARDED one: a tolerance carrying its own resolution floor, and a
        # flat energy read as converged ONLY when the step is still healthy.  Item 6 holds the measurements.
        # AND |grad| IS NOW REPORTED IN HARTREES, which is the only form a reader can act on.  The note above says
        # what the bare number cannot do; the remedy is to convert it, and the conversion needs a curvature.  Along
        # the search direction the energy still reachable is the Newton decrement  dg^2 / (2 <d,Hd>), and <d,Hd>
        # follows from ONE central difference of the exact gradient -- two gradient evaluations, once per run,
        # measured at 2.9 energy evaluations on a correlation layer.  Verified the same day (JAC_EOL_HESSCHECK):
        # the product is symmetric to 9.4e-09 and <d,Hd> agrees with the energy's own second difference to
        # 1.000016.
        #   IT MUST BE THE SEARCH DIRECTION AND NOT THE RAW GRADIENT, which was tried first and is WRONG for a
        # reason worth keeping: the raw gradient lives in the full b-space and therefore carries components along
        # the NEGATIVE-ENERGY branch, where the curvature is enormous and negative.  Measured, stable across four
        # decades of the difference step so it is not round-off: <g,Hg>/<g,g> = -1.53e+04 Ha on the reference layer
        # and -4.03e+03 on the correlation layer, against 2mc^2 = 3.76e+04 Ha.  Used naively that reports "the
        # curvature is not positive, this is not a minimum" at EVERY exit of EVERY run -- a spurious saddle.  The
        # search direction is built by virtualDirections inside the positive branch, and its curvature is positive
        # and physical.  This is also a sharper statement of why |grad| is not a convergence measure: it is not
        # merely unscaled, it is partly supported on rotations the solver is forbidden to make.
        #   IT IS A LOWER BOUND AND IS LABELLED AS ONE: it looks only along one direction, so others may hold more.
        energyStillAvailable = function()
            dHd = 0.
            gg  = 0.;    for sh in activeSubshells   gg = gg + sum( dir[sh].^2 )   end
            gg < 1.0e-30  &&  return( 0.0 )
            epsH = 1.0e-5 / sqrt(gg)                     ## scaled so the probe displacement is ~1e-5 in norm
            bp = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
            bm = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
            for  sh  in  activeSubshells
                bp[sh] = bVectors[sh] + epsH * dir[sh];    bm[sh] = bVectors[sh] - epsH * dir[sh]
            end
            gp = SelfConsistent.computeOrbitalGradient(bp, coeffs1p, coeffs2p, basis.subshells,
                                                               primitives, nucPot, storage, matrixB)
            gm = SelfConsistent.computeOrbitalGradient(bm, coeffs1p, coeffs2p, basis.subshells,
                                                               primitives, nucPot, storage, matrixB)
            for  sh  in  activeSubshells   dHd = dHd + sum( dir[sh] .* (gp[sh] - gm[sh]) ) / (2*epsH)   end
            dHd <= 1.0e-30  &&  return( NaN )            ## no positive curvature along the direction taken
            return( dg*dg / (2*dHd) )
        end
        if  gNorm < settings.accuracyScf
            stopReason = "converged";   println(">> [EOL-C3] CONVERGED at iteration $iter: |grad| = $gNorm < accuracyScf = " *
                    "$(settings.accuracyScf), tStep = $tStep.")
            eAvailAtExit = energyStillAvailable();    tStepAtExit = tStep;    eAtExit = e0
            break
        end

        # AND CONVERGED WHEN THE ENERGY CAN NO LONGER BE RESOLVED -- priority item 6, added 05-Sep-2026.
        # |grad| IS NOT SCALE-FREE: its size is set by the excitation energies that set the curvature, ~1000 Ha
        # at Z = 92 against ~1 Ha at Z = 4, so one threshold cannot serve both.  At Z = 92 the run is DONE and
        # cannot say so -- the energy flat to 4e-15 Ha, a central difference of the functional returning pure
        # round-off, ~3e-09 Ha still available -- while |grad| = 0.0023 reads as unconverged against 1e-6.
        # THIS IS NOT THE NAIVE STATIONARY-ENERGY TEST THAT WAS REMOVED ON 03-Sep, and the difference is the
        # whole design.  That one compared |dE| against accuracyScf, which the energy reaches long before the
        # gradient does, and it stopped every step of Be Scenario A short (iterations 3, 15, 37, 80 against
        # 63, 154, 51, 40).  This one compares |dE| against the RESOLUTION OF THE ENERGY ITSELF -- 32 eps|E|,
        # a few units in the last place -- so it cannot fire while the energy is still measurably falling; on
        # Be that floor is ~1e-13 while |dE| plateaus at 1e-12, and the test stays silent.  It fires only when
        # the arithmetic can no longer tell two successive energies apart.
        # AND IT IS GUARDED BY THE STEP, because a PLATEAU IS NOT A MINIMUM: on Cf^17+ the energy sat
        # stationary to 1e-12 for 55 iterations with tStep at 1e-9, and allowed to run the calculation escaped
        # and fell a further 3.5e-4 Ha, moving the clock transition by 34 cm^-1.  A flat energy measured
        # during a collapsed step says nothing at all.
        eFloor = 32 * eps(abs(e0))
        if  iter > 1  &&  tStep >= stepFloor  &&  abs(e0Prev - e0) < eFloor
            stopReason = "energy below its own resolution"
            println(">> [EOL-C3] CONVERGED at iteration $iter: the energy moved by " *
                    @sprintf("%.2e", abs(e0Prev - e0)) * " Ha, below its own resolution of " *
                    @sprintf("%.2e", eFloor) * " Ha, with a healthy step (tStep = " * @sprintf("%.2e", tStep) *
                    ").  |grad| = $gNorm is reported as a HINT and is not the test: it is not scale-free.")
            eAvail = energyStillAvailable();    eAvailAtExit = eAvail
            tStepAtExit = tStep;    eAtExit = e0
            if  isnan(eAvail)
                println(">> [EOL-C3]   the curvature along the search direction is NOT POSITIVE, so this point is " *
                        "not a minimum along it;  treat the result as a stationary point only.")
            else
                println(">> [EOL-C3]   energy still available along the search direction: " *
                        @sprintf("%.2e", eAvail) * " Ha -- a LOWER BOUND, since only that one direction was " *
                        "measured.  CONVERGED here means stationary to machine resolution FROM THIS START, not " *
                        "the global minimum.")
            end
            break
        end
        e0Prev = e0

        # THE STEP IS INHERITED FROM THE PREVIOUS ITERATION, AND IT MUST BE.  Resetting it to 1.0 here was
        # tried on 31-Aug-2026 and FAILS OUTRIGHT: 24 halvings from unity reach only 1/2^24 = 6.0e-8, while
        # the steps this surface actually accepts at that stage are 4e-9 to 1e-8 -- SMALLER than the search
        # can reach from a unit start -- so no descent is found and the run stops.  Measured on Cf^17+ with
        # an SD layer into {7s,7p}: dead at iteration 11 against 100 iterations of real progress with the
        # step inherited.  The small steps are what the surface requires, not the residue of a collapse.
        # What IS too slow is the recovery rate; see the growth factor at the acceptance below.
        accepted = false
        # ZERO-STEP CONSISTENCY OF THE LINE-SEARCH BASELINE, off unless JAC_EOL_ZEROCHECK is set.  e0 is measured
        # on the projected vectors AFTER restoreFrozen!; a trial is measured on the projected vectors WITHOUT it.
        # If the two disagree at tStep = 0 then `eTrial < e0` compares two different functions and no step can
        # win, however good the direction -- which is the "no descent at iteration 1" seen wherever many
        # subshells are frozen.
        if  haskey(ENV, "JAC_EOL_ZEROCHECK")
            zB = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
            (zProj, _)  = SelfConsistent.projectOntoPositiveBranch(zB, basis.subshells, primitives, nucPot,
                                                        matrixB, storage; spectrum=posSpectrum)
            (_, eZero)  = SelfConsistent.energyFromBVectorsSplit(zProj, coeffs1p, coeffs2p, basis.subshells,
                                                        primitives, grid, nucPot, isFrozenSub, frozenRk)
            zRest       = restoreFrozen!( Dict{Subshell, Vector{Float64}}( sh => copy(zProj[sh])
                                                        for sh in basis.subshells ) )
            (_, eZeroR) = SelfConsistent.energyFromBVectorsSplit(zRest, coeffs1p, coeffs2p, basis.subshells,
                                                        primitives, grid, nucPot, isFrozenSub, frozenRk)
            @printf(">> [EOL-ZERO] iter %d: e0 = %.12f  eTrial(0) = %.12f  diff = %+.3e  restored-diff = %+.3e  nFrozen = %d\n",
                    iter, e0, eZero, eZero - e0, eZeroR - e0, length(frozenSubshells))
            flush(stdout)
        end
        # THE STEP SCAN, 07-Sep-2026 -- this is what found the sign discontinuity, so it is kept.
        # It walks the trial step over four orders around the current one and reports the energy together with
        # the cosine of each TABULATED orbital against its value at the smallest step.  A DISCONTINUOUS
        # functional shows here and almost nowhere else: on C-like U it gave a jump of exactly 0.032 Ha between
        # 2.20e-08 and 2.25e-08, with one orbital's cosine going to -1 while every other diagnostic in this
        # file -- the b-vectors, negW, the orthonormality deviation, mtp -- was unchanged to five digits.
        # COMPARE THE TABULATED ORBITAL AND NOT THE VECTOR: the sign convention is applied inside
        # generateOrbitalFromVector, so a flip is invisible in the b-vector and invisible in wSign itself.
        # Off unless JAC_EOL_BISECT names an iteration.
        if  get(ENV, "JAC_EOL_BISECT", "") == string(iter)
            @printf(">> [EOL-BISECT] iteration %d, e0 = %.12f\n", iter, e0);   flush(stdout)
            bisectRef = Dict{Subshell, Vector{Float64}}()
            for  ts  in  tStep .* [0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0, 4.0]
                bB = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in basis.subshells )
                for  sh  in  activeSubshells
                    v = bVectors[sh] + ts * dir[sh]
                    bB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
                end
                # the SAME functional on the UNPROJECTED vectors, to place the jump on one side of the
                # projection or the other
                bRaw = Dict{Subshell, Vector{Float64}}( sh => copy(bB[sh])  for sh in basis.subshells )
                restoreFrozen!(bRaw)
                (_, eRaw) = SelfConsistent.energyFromBVectorsSplit(bRaw, coeffs1p, coeffs2p, basis.subshells,
                                                        primitives, grid, nucPot, isFrozenSub, frozenRk)
                (bProj, bNeg) = SelfConsistent.projectOntoPositiveBranch(bB, basis.subshells, primitives,
                                                        nucPot, matrixB, storage; spectrum=posSpectrum)
                restoreFrozen!(bProj)
                (_, bE) = SelfConsistent.energyFromBVectorsSplit(bProj, coeffs1p, coeffs2p, basis.subshells,
                                                        primitives, grid, nucPot, isFrozenSub, frozenRk)
                # the TABULATED orbital, not the b-vector: generateOrbitalFromVector canonicalises the sign
                # with sum(P[1:30]), which can cross zero and flip the whole orbital -- invisible in the vector
                # compare the TABULATED orbital against the one at the smallest step: a sign flip inside
                # generateOrbitalFromVector is invisible in wSign, because it is applied before the return
                ovl = join([ begin
                        ob = Bsplines.generateOrbitalFromVector(sh, 0.0, bProj[sh], primitives)
                        rf = get!(bisectRef, sh, ob.P)
                        n  = min(length(ob.P), length(rf))
                        d  = sum(ob.P[1:n] .* rf[1:n]) / sqrt(sum(rf[1:n].^2) * sum(ob.P[1:n].^2))
                        @sprintf("%s %+.4f", string(sh), d)
                    end  for sh in basis.subshells ], " ")
                @printf(">> [EOL-BISECT] tStep %.4e  dE %+.6e  | cos to reference: %s\n", ts, bE - e0, ovl)
                flush(stdout)
            end
        end
        for  trial = 1:24
            newB = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in basis.subshells )
            for  sh  in  activeSubshells
                v = bVectors[sh] + tStep * dir[sh]
                newB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
            end
            # The projection must happen BEFORE the acceptance test, not after it.  Projecting and
            # re-orthonormalizing moves the orbitals, so accepting `newB` on the strength of its own energy
            # and then storing the PROJECTED vector stores something that was never tested -- and the
            # projection gives a little of the gain back each step.  That broke monotonicity: on Li the
            # driver reduced the gradient 19-fold while the energy ROSE by 3.7e-7 Ha, which a descent
            # method cannot do.  Testing the projected vector restores  E(new) < E(old)  by construction,
            # and with it the guarantee that the CI eigenvalue falls too (it is bounded above by this
            # fixed-coefficient functional, and equals it at the previous orbitals).
            # THE PROJECTION STAYS INSIDE THE ACCEPTANCE TEST, and letting the search leave the branch was
            # MEASURED AND REJECTED on 02-Sep-2026, after the frozen-orbital repair of the same day removed the
            # earlier objection to that measurement.  Judging a trial on its UNPROJECTED energy and storing the
            # projected vector banks a real loss for a fake gain, and the two compound: on Be-like Z = 92 the
            # reference layer went from -12040.8576547019 to -11691.0209575892, i.e. 350 Ha UPHILL, with |grad|
            # exploding from 1.52 to 11195.9 by the second iteration.  That is variational collapse into the
            # negative-energy sea, and it is what this projection exists to prevent.
            (projB, negW) = SelfConsistent.projectOntoPositiveBranch(newB, basis.subshells,
                                                    primitives, nucPot, matrixB, storage; spectrum=posSpectrum)
            restoreFrozen!(projB)
            (_, eTrial) = SelfConsistent.energyFromBVectorsSplit(projB, coeffs1p, coeffs2p, basis.subshells,
                                                 primitives, grid, nucPot, isFrozenSub, frozenRk)
            # THE SHAPE OF THE ENERGY ALONG THE SEARCH DIRECTION, off unless JAC_EOL_LINEPROFILE is set.  When a
            # line search fails, the three possible causes are told apart by this one picture: an energy that
            # RISES at every trial means the direction is uphill (the gradient is wrong); one that stays FLAT at
            # the 1e-13 level means the step is not moving the orbitals at all (the projection is giving it
            # back); one that FALLS and is still rejected means the acceptance test is misaccounting.  `moved`
            # is the size of the displacement that actually survived projection and renormalization.
            if  haskey(ENV, "JAC_EOL_LINEPROFILE")
                moved = 0.
                for  sh  in  activeSubshells   moved = moved + sum( (projB[sh] - bVectors[sh]).^2 )   end
                # <grad, realized displacement> against tStep*<grad,dir>: the derivative along the step the
                # search ACTUALLY takes, versus the one it planned.  They differ by whatever the renormalization
                # and the positive-branch projection do to the direction.
                dgReal = 0.
                for  sh  in  activeSubshells   dgReal = dgReal + sum( grad[sh] .* (projB[sh] - bVectors[sh]) )   end
                @printf(">> [EOL-LINE] iter %2d trial %2d: tStep = %.3e  eTrial-e0 = %+.6e  moved = %.3e  planned = %+.4e  realized = %+.4e  negW = %.3e\n",
                        iter, trial, tStep, eTrial - e0, sqrt(moved), tStep*dg, dgReal, negW)
                flush(stdout)
            end
            if  eTrial < e0
                bVectors = projB
                if  printout  &&  negW > 1.0e-8
                    println(">> [EOL-C3] removed negative-branch weight $negW from the step.")
                end
                # THE 1.3x IS SLOW ON PURPOSE-BY-ACCIDENT, AND IT MUST STAY UNTIL THE STATIONARY-ENERGY TEST
                # IS SOUND.  Growing 4x on a first-trial success was tried on 31-Aug-2026 and made the result
                # WORSE, for a reason worth recording.  A collapse takes log(1e7)/log(1.3) = 61 iterations to
                # climb back, and on Cf^17+ (SD layer into {7s,7p}) that showed as a 55-iteration plateau with
                # |grad| frozen at 0.007384 -- so growing faster looks obviously right.  But during that
                # plateau the ENERGY is stationary to 1e-12 while the calculation is NOT converged: allowed to
                # run, it escapes and falls a further 3.5e-4 Ha, moving the clock transition by 34 cm^-1.  The
                # stationary-energy exit below is blocked during the plateau only because tStep sits at 1e-9,
                # BELOW its stepFloor guard.  Growing 4x lifts the step over that floor while the energy is
                # still flat, the exit fires at iteration 17, and the run returns the less converged answer
                # (8834.25 cm^-1 against 8867.97 after 100 iterations).  So the slow growth is compensating
                # for an exit test that cannot tell a plateau from a minimum.  FIX THE TEST FIRST; the growth
                # factor is then free to be raised, and should be.
                # THE ARMIJO RATIO: what the step actually bought, over what the direction promised.
                # predicted decrease = -tStep * <grad, dir> (dg < 0 by the descent guard above), actual =
                # e0 - eTrial.  A ratio near 1 means the linear model describes this surface and the small
                # steps are the surface's own curvature; a ratio near 0 means the direction is poor and the
                # halving is the search compensating for it.  The two call for opposite repairs, and nothing
                # in the solver reported this quantity before, which is why three attempts to cure the
                # plateau of 30/31-Aug-2026 were aimed at the symptom.
                # TWO ratios, because the line search's model and its measurement are taken at DIFFERENT
                # points.  ratio1 uses tStep*<grad,dir>, the decrease the direction promises for the step as
                # PLANNED.  But the point actually evaluated is not b + tStep*dir: each orbital is
                # renormalized and then projected onto the positive branch first.  ratio2 therefore uses the
                # displacement that really happened, <grad, projB - b>.  For a genuine gradient ratio2 must
                # tend to 1 as the step shrinks -- that is what a derivative means -- so if ratio1 sits at a
                # constant far from 1 while ratio2 approaches it, the model is describing a step the search
                # does not take, and the repair is to model the displacement instead.  If BOTH stray, the
                # gradient itself disagrees with the functional and the fault is in computeOrbitalGradient.
                if  printout
                    predicted1 = -tStep * dg
                    predicted2 = 0.
                    for  sh  in  activeSubshells
                        predicted2 = predicted2 - sum( grad[sh] .* (projB[sh] - bPrev[sh]) )
                    end
                    r1 = predicted1 > 0. ? (e0 - eTrial)/predicted1 : NaN
                    r2 = predicted2 > 0. ? (e0 - eTrial)/predicted2 : NaN
                    @printf(">> [EOL-C3]    accepted at trial %2d, tStep = %.3e, planned = %+.4f, actual-disp = %+.4f\n",
                            trial, tStep, r1, r2)
                end
                # THE GROWTH FACTOR STAYS AT 1.3, AND THAT WAS MEASURED RATHER THAN ASSUMED.  It only matters
                # after a collapse -- recovery from 1e-8 takes log(1e7)/log(1.3) = 61 iterations -- which was a
                # real cost while the directional derivative was five to nine times too steep and collapses
                # were routine (items 121, 122).  With the derivative exact the steps run at 0.03 to 1.0 and
                # rarely collapse, so the rate stopped being the binding constraint.  Measured 31-Aug-2026 on
                # Be Scenarios A and B, growth applied only on a first-trial success:
                #     1.3  all six A steps converge at 63, 154, 51, 40, 51, 40
                #     2.0  all six converge at 39, 145, 54, 41, 54, 41 -- faster on two, slower on two
                #     4.0  A REGRESSION: default-grid step 1 stops on a stagnant gradient at 57 where 1.3
                #          converges at 63; B stops stagnant at 46
                # 2.0 is a wash and 4.0 loses a convergence, so there is no case for changing a working
                # default.  A faster step that overshoots into "no descent" trades a real convergence for a
                # marginal speedup, which is the wrong bargain in a solver whose problem was trustworthiness.
                accepted = true;    tStep = min(1.0, 1.3*tStep);    break
            end
            tStep = tStep / 2
        end
        if  printout
        end
        if  !accepted  &&  method == :lbfgs  &&  !isempty(sHist)
            # A line search that finds no descent along an L-BFGS direction means the stored curvature is
            # no longer describing this surface -- unsurprising, since virtualDirections is rebuilt every
            # iteration and the pairs then mix vectors from different subspaces. Discard the history and
            # retry from the preconditioned gradient rather than giving up.
            if  printout    println(">> [EOL-C3] L-BFGS history discarded at iteration $iter; retrying.")   end
            empty!(sHist);   empty!(yHist);   empty!(rhoHist);   tStep = 1.0
            for  trial = 1:24
                newB = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in basis.subshells )
                for  sh  in  activeSubshells
                    v = bVectors[sh] + tStep * sVec[sh]
                    newB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
                end
                (projB, negW) = SelfConsistent.projectOntoPositiveBranch(newB, basis.subshells,
                                                    primitives, nucPot, matrixB, storage; spectrum=posSpectrum)
                restoreFrozen!(projB)                    ## measure the point that would be stored; see the main loop
                (_, eTrial) = SelfConsistent.energyFromBVectorsSplit(projB, coeffs1p, coeffs2p, basis.subshells,
                                                         primitives, grid, nucPot, isFrozenSub, frozenRk)
                if  eTrial < e0    bVectors = projB;   accepted = true;   break    end
                tStep = tStep / 2
            end
        end
        if  !accepted
            eAvailAtExit = energyStillAvailable();    tStepAtExit = tStep;    eAtExit = e0
            stopReason = "no descent";   println(">> [EOL-C3] STOPPED at iteration $iter: no descent found along the search direction " *
                    "(tStep fell to $tStep), with |grad| = $gNorm.  This is NOT convergence.")
            Defaults.warn(AddWarning(), "SelfConsistent.solveOptimizedLevelFieldByRotation(): the EOL field did NOT " *
                          "converge -- no descent found at iteration $iter, |grad| = " * @sprintf("%.1e", gNorm) *
                          ".  The energies are NOT self-consistent.")
            break
        end
        # A STAGNANT ENERGY IS ONLY CONVERGENCE IF THE STEP IS STILL HEALTHY.  The test compares successive
        # CI eigenvalues, and a variational energy is STATIONARY at a minimum: |dE| falls as the SQUARE of the
        # orbital error while |grad| falls linearly, so |dE| < 1e-11 is reached long before convergence
        # whenever the steps have become small -- and then it reports a collapsed line search as a converged
        # calculation.  Measured at the moment it used to fire: tStep = 4.1e-9 on a Be RAS correlation layer
        # and 3.5e-11 on a carbon one, i.e. down eight to eleven orders from unity, with |grad| still 0.037
        # and 0.998.  The collapse is TEMPORARY -- allowed to continue, carbon recovers a step of 0.125 and
        # converges at |grad| = 4.7e-6, sixty-five mHa BELOW where it used to stop.  So the remedy is not a
        # smaller threshold but the extra condition: stagnation ends the iteration only when the step that
        # produced it was of usable size.
        # THE ITERATION ENDS WHEN THE GRADIENT STOPS IMPROVING -- NOT WHEN THE ENERGY STOPS MOVING.
        # A stationary energy is a weak signal near a minimum: the energy is quadratic there and flattens long
        # before the gradient does, which is why |dE| reaches 1e-12 while |grad| is still 1e-5.  Ending on it
        # stopped every step of Be Scenario A short of convergence -- at iterations 3, 15, 37 and 80 -- while
        # the SAME runs with the exit disabled converged at 63, 154, 51 and 40 with |grad| < 1e-6.  No
        # threshold on |grad| separates those cases (6.95e-5 still converges if allowed), so the test itself
        # was the wrong one.
        # The question this exit exists to answer is whether more iterating will help, and the quantity that
        # answers it is the gradient's own progress.  It is given `stagnationWindow` iterations to beat its
        # best value by any margin; if it cannot, no further iterating will, and the run stops saying so.
        # The original justification for the energy test -- that |grad| plateaus at a floor set by the basis
        # and the projection -- was measured while the directional derivative was five to nine times too steep
        # (items 121 and 122).  With the gradient exact there is no such floor.
        # THE BUDGET CASE IS MEASURED HERE, on the last iteration, because `dir` and the closure that needs it are
        # local to this loop and are gone by the time the budget is found to be exhausted.  One extra measurement
        # on one iteration, and it is the case a user most wants a number for: a run that simply ran out.
        if  iter == Basics.maxIterations(settings.scfRoute)
            eAvailAtExit = energyStillAvailable();    tStepAtExit = tStep;    eAtExit = e0
        end
        # THE INCREMENT HAS STOPPED MOVING UNDER A DOUBLING OF THE EFFORT -- the one exit that is a convergence
        # criterion for a CORRELATION LAYER, off unless the route asks for it (incrementTolerance > 0).
        #
        # WHY NOT THE ENERGY CHANGE PER ITERATION, which is the obvious thing and is wrong.  A correlation layer
        # creeps: measured on Be Scenario A, a stationary-energy test stopped it at iterations 3, 15, 37 and 80
        # where the same runs were still improving at 63, 154, 51 and 40.  A one-iteration difference cannot tell
        # a layer that has converged from one that is descending slowly, because both look flat from step to step.
        #
        # WHAT WORKS IS A LONG BASELINE.  Compare the descent achieved so far, D(n) = E(1) - E(n), against the
        # descent at HALF the iterations, D(n/2).  If doubling the effort from n/2 to n changed the layer's
        # increment by less than the tolerance, doubling it again will not change it either, and the layer is as
        # converged as it is going to get.  This is the doubling test of the RAS strategy, applied continuously
        # instead of by running the whole computation twice -- the maintainer's suggestion, 07-Oct-2026 -- so the
        # budget becomes a CEILING rather than the thing that decides the answer.
        #
        # THE SCALE IS THE DESCENT ITSELF, which is what makes the test free of units and of Z.  An absolute
        # tolerance cannot serve both a layer worth 160 mHa (Ti III 3d^2) and one worth 0.06 mHa (Ca+ 4s with core
        # s singles), and those two were measured on the same afternoon.
        #
        # NOT BELOW SIX ITERATIONS, measured: at a budget of 3 a reference step returned an energy BELOW its own
        # converged value, so a test allowed to fire that early can stop on noise.
        if  settings.scfRoute isa Basics.RotationRoute  &&  Basics.incrementTolerance(settings.scfRoute) > 0.  &&
                    iter >= 6  &&  length(eHist) == iter
            half = div(iter, 2)
            dNow = eHist[1] - eHist[iter];      dHalf = eHist[1] - eHist[half]
            if  dNow > 0.  &&  abs(dNow - dHalf) <= Basics.incrementTolerance(settings.scfRoute) * dNow
                eAvailAtExit = energyStillAvailable();    tStepAtExit = tStep;    eAtExit = e0
                stopReason = "increment stable under doubling"
                iterDone   = iter
                println(">> [EOL-C3] CONVERGED at iteration $iter: the layer's increment is " *
                        @sprintf("%.6e Ha", dNow) * " and was " * @sprintf("%.6e Ha", dHalf) *
                        " at iteration $half, so DOUBLING the effort moved it by " *
                        @sprintf("%.2f %%", 100*abs(dNow-dHalf)/dNow) * ", within the requested " *
                        @sprintf("%.2f %%", 100*Basics.incrementTolerance(settings.scfRoute)) *
                        ".  |grad| = $gNorm, tStep = $tStep.")
                break
            end
        end
        stagnationWindow = 20
        # THE ITERATION NOW ENDS WHEN THE ENERGY STOPS IMPROVING, NOT WHEN THE GRADIENT DOES -- 07-Sep-2026.
        # This reverses the decision of 03-Sep, and it reverses it because the surface changed underneath it.
        # THEN: |grad| was pinned by a discontinuity (the orbital sign could flip mid-line-search, 559c3ea), so
        # a stationary energy really did mean "stuck" and the gradient was the better signal.
        # NOW: the functional is differentiable and the analytic gradient matches a central difference to
        # 1.0000, and on that surface the ordering is the other way round.  MEASURED on C-like U layer 2 with
        # full steps: |grad| falls to 7.77e-02 by iteration 26 and RISES to 1.57e-01 by 46, while the ENERGY
        # falls steadily by ~0.025 Ha over the same span -- so a gradient-watching exit stops a run that is
        # descending.  That is exactly item 6's Z = 4 symptom, "both exits fire while the energy is STILL
        # FALLING", and it is what the item prescribed fixing: progress on the TOTAL ENERGY as GRASP does,
        # with the gradient reported as a HINT.
        # AND THE ENERGY WATCHED IS THE CI EIGENVALUE, NOT THE FUNCTIONAL e0.  They are not interchangeable and
        # the difference decides runs: the EOL functional is REBUILT every iteration, because its angular
        # coefficients come from the current CI vector, so comparing e0 across iterations compares DIFFERENT
        # FUNCTIONS.  Measured on C-like U layer 2: e0 sat at -1.1418373642578115 from iteration 29 while the
        # CI energy fell from -14351.7287 to -14351.7457, i.e. 0.017 Ha of real descent that an e0-watching
        # exit reads as stagnation.  The CI eigenvalue is the quantity that is the SAME question every
        # iteration, and it is the one a user is given.
        # THE MARGIN IS THAT ENERGY'S OWN RESOLUTION, not a tolerance: an improvement smaller than the
        # arithmetic can represent is not an improvement.  It is NOT capped here as it is in the convergence
        # exit -- there a plateau must not be read as convergence, whereas here reading it as stagnation is
        # exactly right.
        eCI = multiplet.levels[1].energy
        if  eCI < bestE - 32 * eps(abs(eCI))    bestE = eCI;    bestEIter = iter    end
        if  gNorm < bestGNorm * (1.0 - 1.0e-3)    bestGNorm = gNorm;    bestGIter = iter    end
        # THE GUARD THE NOTE ABOVE PRESCRIBES, WIRED IN 05-Sep-2026 (priority item 6).  `stepFloor` was
        # defined here and never used: the comment said "stagnation ends the iteration only when the step that
        # produced it was of usable size" and the condition did not test it, so a TEMPORARY step collapse --
        # the carbon case, 3.5e-11 with |grad| still 0.998 -- could end the run 65 mHa short of where it goes
        # when allowed to continue.  The guard can only let a run continue, never stop it earlier.
        # A COLLAPSE THAT NEVER RECOVERS IS NOT A TEMPORARY ONE -- 06-Sep-2026, and this bounds the guard added
        # in d0d66fc.  That guard exists because a collapse CAN be temporary: a carbon case fell to 3.5e-11 with
        # |grad| still 0.998 and, allowed to continue, recovered a step of 0.125 and converged 65 mHa lower.  But
        # unbounded it turns a dead run into a full iteration budget: MEASURED on C-like U at 1e-8, the step
        # collapsed to ~2e-08 at iteration 47 and stayed there for 350 more iterations with |grad| frozen at
        # 0.0770 and the energy flat in its eleventh decimal, while the guard blocked the exit at every one.
        # So the step is given `stagnationWindow` iterations to recover, the same budget the gradient gets.
        if  tStep >= stepFloor    collapsedSince = 0
        elseif  collapsedSince == 0    collapsedSince = iter
        end
        stepIsDead = collapsedSince > 0  &&  iter - collapsedSince >= stagnationWindow
        if  haskey(ENV, "JAC_EOL_STAGDIAG")  &&  iter - bestGIter >= stagnationWindow
            @printf(">> [EOL-STAG] iter %3d  bestGIter %3d  tStep %.3e  collapsedSince %3d  -> %s\n",
                    iter, bestGIter, tStep, collapsedSince,
                    (tStep >= stepFloor || stepIsDead) ? "EXIT" : "held open");  flush(stdout)
        end
        if  iter - bestEIter >= stagnationWindow  &&  (tStep >= stepFloor || stepIsDead)  &&
                                                      !haskey(ENV, "JAC_EOL_NOSTATEXIT")
            eAvailAtExit = energyStillAvailable();    tStepAtExit = tStep;    eAtExit = e0
            stopReason = "energy stagnated";   println(">> [EOL-C3] stopped at iteration $iter: the ENERGY has " *
                    "not improved on $bestE since iteration $bestEIter, $stagnationWindow iterations ago" *
                    (stepIsDead ? ", and the step has been collapsed below $stepFloor since iteration $collapsedSince" : "") * ".  " *
                    "A converging UPPER BOUND.  |grad| = $gNorm is a HINT and not the test: it is not scale-free, " *
                    "and on a smooth surface it can RISE while the energy falls (best was $bestGNorm at $bestGIter).")
            Defaults.warn(AddWarning(), "SelfConsistent.solveOptimizedLevelFieldByRotation(): the EOL field stopped " *
                          "on a stagnant ENERGY at iteration $iter with |grad| = " * @sprintf("%.1e", gNorm) *
                          ".  The energy is a converging UPPER BOUND.")
            break
        end
    end
    if  stopReason == ""
        println(">> [EOL-C3] STOPPED after $iterDone iterations: the budget of $(settings.scfRoute) was reached " *
                "with |grad| = $gNorm and tStep = $tStep.  This is NOT convergence: the rotation route descends " *
                "monotonically, so a larger maxIterations is very likely to keep improving the energy.")
        Defaults.warn(AddWarning(), "SelfConsistent.solveOptimizedLevelFieldByRotation(): the EOL field did NOT " *
                      "converge -- the budget of $(Basics.maxIterations(settings.scfRoute)) iterations was reached with " *
                      "|grad| = " * @sprintf("%.1e", gNorm) * ".  The energies are NOT self-consistent.")
    end

    # THE VERDICT IS RECORDED SO A CALLER CAN ASK INSTEAD OF GREPPING, which is what priority item 40 asked for and
    # all it asked for: the checks already existed and already printed, and none of them returned anything a script
    # could assert on.  The return value of this function is UNCHANGED -- it is still the Multiplet -- so nothing
    # downstream has to be touched; the verdict is fetched with SelfConsistent.lastScfVerdict().
    #   ADVISORY, NEVER GATING.  Nothing in JAC reads this, and no computation is blocked by it.  A surprising
    # physical result must remain obtainable with every field of it looking wrong, which is the maintainer's
    # constraint on the whole item and the reason it is a record rather than a check.
    descentNow  = length(eHist) >= 2         ? eHist[1] - eHist[end]            : 0.
    descentHalf = length(eHist) >= 4         ? eHist[1] - eHist[div(length(eHist),2)] : 0.
    GBL_EOL_LAST_VERDICT[] = ScfVerdict(stopReason == "converged" || stopReason == "energy below its own resolution" ||
                                        stopReason == "increment stable under doubling",
                                        stopReason, iterDone, eAvailAtExit, gNorm, tStepAtExit, eAtExit,
                                        descentNow, descentHalf)

    # THE REQUESTED INTERACTION IS APPLIED HERE, ONCE, ON THE CONVERGED ORBITALS.  Until 01-Sep-2026 this
    # function returned the multiplet built inside its own iteration by the EOL machinery
    # (cacheCsfPairCoefficientsEOL / buildCIMatrixEOL), which is PURE COULOMB BY CONSTRUCTION -- so
    # settings.eeInteractionCI was silently ignored on the whole EOL path, and a Breit or QED request
    # returned Coulomb numbers with no warning.  Its sibling solveOptimizedLevelField always ended with this
    # call; performSCF dispatches here, to the one that dropped it, and the file's own note above (that the
    # non-Coulomb terms "are added only once, at the final Hamiltonian.performCIKinkAware call") described
    # the sibling and not this function.
    #   THE SCF LOOP ITSELF STAYS PURE COULOMB, deliberately: that is the variational functional the rotation
    # minimizes, and adding Breit inside the iteration would change what is being optimized rather than what
    # is being reported.  The correction belongs at the end, which is also where AL and DFS apply it.
    # EVERY ORBITAL IS GIVEN ITS DIAGONAL LAGRANGE MULTIPLIER eps_a = <a|F_a|a>, 12-Sep-2026.  Until then the
    # energy field was 0.0 here, which was HONEST rather than an oversight -- a rotation-optimised orbital is
    # not the eigenfunction of any one-particle operator, so it has no eigenvalue to carry -- but it cost
    # something real: omega = factor |E_a - E_c| / c came out ZERO for every pair of such orbitals, so a
    # frequency-dependent CoulombBreit(factor > 0) silently returned its own omega -> 0 limit on any EOL or
    # RAS basis, and anything else reading an orbital energy was equally stuck.
    #   THE MULTIPLIER IS THE GENERALISATION THE ORBITAL DOES HAVE.  Where the orbital IS an eigenvector of
    # F_a -- which is exactly the AL case -- the Rayleigh quotient and the stored eigenvalue are the same
    # number, so this defines the same quantity AL already reports rather than a second convention.
    #   The one-particle expectation <a|h_D|a> is NOT this quantity and must not be substituted: it omits the
    # electron-electron interaction and comes out 1.7x to 12x too deep while looking well-behaved (measured,
    # tools/probe-eolOrbitalEnergy.jl).
    #   The cost is one Fock build per subshell, once, after convergence.
    #   THE TARGET LEVELS ARE RE-SELECTED ON THE FINAL ORBITALS, not carried out of the loop: the multiplier is a
    # property of the converged orbital set, and the loop's last iterate was built on the previous one.
    finalTmpOrbs = Dict{Subshell, Orbital}()
    for  sh  in  basis.subshells
        finalTmpOrbs[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives; canonicalize=false)
    end
    finalTmpBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, finalTmpOrbs)
    finalLevels   = Level[]
    let  radial1pF = Dict{Tuple{Subshell,Subshell},Float64}(),
         radial2pF = Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}(),
         vkF       = Dict{Tuple{Int64,Subshell,Subshell,Int64},Vector{Float64}}()
        for  sym  in  relevantSyms
            cache = blockCaches[sym];    idxCsf = cache.idxCsf
            mtx = SelfConsistent.buildCIMatrixEOL(cache, finalTmpOrbs, grid, nucPot,
                                                  radial1pF, radial2pF, vkF)
            append!( finalLevels, SelfConsistent.diagonalizeBlockEOL(sym, idxCsf, mtx, finalTmpBasis) )
        end
    end
    finalTargetLevels = SelfConsistent.selectTargetLevelsEOL(Basics.sortByEnergy(Multiplet("EOL-ByRotation",
                                                             finalLevels)), settings.levelSelectionCI)
    # The off-diagonal CSF-pair terms follow the SAME convention the sibling solver refines with, so that the
    # two report the same quantity; see the note at `unscaledOff` there.
    unscaledOffF = SelfConsistent.GBL_EOL_UNSCALED_OFFDIAGONAL
    if  unscaledOffF
        (_, coeffs2pF)    = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, finalTargetLevels; pairs=:diagonal)
        (_, coeffs2pFOff) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, finalTargetLevels; pairs=:offdiagonal)
    else
        (_, coeffs2pF)    = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, finalTargetLevels)
        coeffs2pFOff      = Coefficient2p[]
    end
    genOccF     = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, finalTargetLevels, basis)
    orbEnergies = SelfConsistent.computeOrbitalEnergiesEOL(basis.subshells, bVectors, coeffs2pF, genOccF,
                                                           primitives, grid, nucPot, storage, matrixB;
                                                           coeffs2pUnscaled=coeffs2pFOff)
    finalOrbitals = Dict{Subshell, Orbital}()
    for  sh  in  basis.subshells
        finalOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, orbEnergies[sh], bVectors[sh], primitives)
    end
    finalBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, finalOrbitals)
    multiplet  = Hamiltonian.performCIKinkAware(finalBasis, nuclearModel, grid, settings; printout=printout)

    return( multiplet )
end


"""
`SelfConsistent.solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                         settings::AsfSettings; printout::Bool=true)`
    ... solves the self-consistent field for the extended-optimal-level (EOL) functional: orbitals are
        optimized against a statistically-(2J+1)-weighted combination of one or more target ASF levels,
        with the combination weights coming from the levels' own CI mixing coefficients -- never a
        user-supplied weight. The target level(s) are selected EXCLUSIVELY via settings.levelSelectionCI
        (either indices or symmetries, never both; if inactive/empty, defaults to indices=[1], a genuine
        OL/single-level computation). This nests a CI diagonalization inside the AL SCF loop
        (GRASP rmcdhf90's scf.f90, algorithm 5.1 of Froese Fischer, Comput. Phys. Rep. 3 (1986) 290):
        diagonalize -> build generalized coefficients from the target levels' mixing vectors -> refine
        every orbital (reusing computeFockMatrix/computeTwoElectronV unchanged) ->
        re-diagonalize -> repeat, converging on both orbital self-consistency and the weighted-average
        energy. Per-CSF-pair angular coefficients (cacheCsfPairCoefficientsEOL) are computed once and
        reused every outer iteration, since they depend only on the fixed CSF list, never on the current
        orbitals or mixing coefficients. A (new) multiplet::Multiplet is returned.

        KNOWN LIMITATION (confirmed 28-Jul-2026, NOT yet fixed): when two or more CSFs of the SAME symmetry
        block compete for the same correlation channel (e.g. Be's 1s^2 2p_1/2^2 vs 1s^2 2p_3/2^2, both
        correlating with 1s^2 2s^2), this implementation can converge to a spurious, winner-take-all fixed
        point where one competing CSF's mixing coefficient is driven to ~0 while a comparably-important
        partner is not -- confirmed against a DFS-Field reference occupying the same CSF space, which lands
        both more bound AND with both CSFs contributing substantially. Root cause: the off-diagonal
        (CSF-pair) coefficients folded into computeTwoElectronV's two-electron potential scale
        LINEARLY in a shrinking CSF's own mixing coefficient, while computeGeneralizedOccupationEOL's
        occupation (the (1.0/occ) divisor in computeFockMatrix) scales QUADRATICALLY in it -- so the
        ratio diverges as that coefficient shrinks, rather than settling. This is the concrete manifestation
        of the "DA/inhomogeneous-term mechanism" gap vs. GRASP's setcof.f90 (which treats within-level
        off-diagonal coupling as a separate inhomogeneous/source term, not folded into the same per-orbital
        homogeneous eigenvalue division) -- see project_eol_implementation.md. Flooring `occ` before the
        division was tried and REJECTED as a fix (non-monotonic in the floor constant).

        MEASURED 09-Sep-2026, and the diagnosis above is CONFIRMED but is only HALF of the error. Be
        1s^2 2s^2 + 1s^2 2p^2 at Z = 4, every route started from the same average-level basis:

            Fock, off-diagonal terms SCALED by 1/occ   E = -14.594824221   2p_3/2 <r> = 10.08, weight 1.9e-05
            Fock, off-diagonal terms UNSCALED          E = -14.610656334   2p_3/2 <r> =  3.41, weight 0.053
            average level                              E = -14.613805206
            rotation route                             E = -14.619514867   2p_1/2, 2p_3/2 agree to 5 digits

        Keeping the off-diagonal contributions OUT of the 1/occ scaling -- Basics.FockRoute's
        unscaledOffDiagonal, which is its default -- removes the winner-take-all collapse and is worth 15.8 mHa.

        WHAT REMAINS IS A DEFECT OF THE EQUATIONS, not of the search, and this is the measurement that shows it:
        STARTED AT the rotation minimum, this solver walks UPHILL by 8.86 mHa to the very same fixed point it
        reaches from an average-level start (-14.610656351 against -14.610656334). Its fixed point is therefore
        not a stationary point of the EOL functional. The reason is that an off-diagonal element such as
        R^k(2s,2s,2p,2p) differentiates, with respect to the 2p orbital, into a term proportional to P(2s) -- an
        inhomogeneous SOURCE on the right-hand side -- while computeTwoElectronV applies it as a matrix
        multiplying P(2p). Removing the 1/occ scaling put that term at the right SIZE; it is still on the wrong
        SIDE, which is what GRASP's setcof.f90 carries as a source term. Until that exists this route converges
        quickly to an answer reliably above the rotation route's, and the signature to watch for is the two
        spin-orbit partners of a subshell disagreeing in mean radius: 2.95 against 3.41 here, where the rotation
        route gives 2.502 and 2.502.
"""
function solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                  settings::AsfSettings; printout::Bool=true)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid

    if  !settings.levelSelectionCI.active  ||  ( isempty(settings.levelSelectionCI.indices) && isempty(settings.levelSelectionCI.symmetries) )
        println(">> [EOL] No levelSelectionCI given; defaulting to indices = [1] -- a genuine OL (single lowest level) computation.")
    elseif  !isempty(settings.levelSelectionCI.indices)  &&  !isempty(settings.levelSelectionCI.symmetries)
        error("stop a; settings.levelSelectionCI must specify EITHER indices OR symmetries for the EOL scheme, not both.")
    end

    # Determine which symmetry block(s) need to be (re-) diagonalized every outer iteration: an explicit
    # symmetries-only selection only ever needs those blocks; index-based (or default) selection needs
    # every block present in the basis to resolve the global, energy-sorted ordering.
    if  settings.levelSelectionCI.active  &&  !isempty(settings.levelSelectionCI.symmetries)
        relevantSyms = unique( settings.levelSelectionCI.symmetries )
    else
        relevantSyms = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
    end

    # (1) Initialize storage and important arrays; determine nuclear potential and mean occupation once --
    # identical to solveAverageLevelField
    if  printout    println(">> [EOL] (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    nucPot  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)

    bVectors = Dict{Subshell, Vector{Float64}}()
    for  sh  in  basis.subshells
        bVectors[sh] = Bsplines.fitVectorToPrimitives(basis.orbitals[sh], primitives, matrixB)
    end
    orbitals = basis.orbitals

    # (2) Cache the (orbital-independent) per-CSF-pair angular coefficients once for every relevant block
    if  printout    println(">> [EOL] Caching per-CSF-pair angular coefficients for symmetries $(relevantSyms) ...")    end
    blockCaches = Dict{LevelSymmetry, SelfConsistent.PairCoefficientCache}()
    for  sym  in  relevantSyms
        blockCaches[sym] = SelfConsistent.cacheCsfPairCoefficientsEOL(sym, basis)
    end

    # (3) Initial CI diagonalization (starting orbitals) to get the first mixing vectors for the target level(s)
    function diagonalizeAllBlocks(currentOrbitals::Dict{Subshell, Orbital})
        tempBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, currentOrbitals)
        # One shared radial-integral cache per outer iteration, reused across every block: the SAME
        # subshell-labeled integral (e.g. "1s-1s") often recurs across many CSF pairs and even across
        # different symmetry blocks, but depends only on currentOrbitals, not on which block/pair asked for it.
        radial1pCache = Dict{Tuple{Subshell,Subshell},Float64}()
        radial2pCache = Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}()
        vkCache       = Dict{Tuple{Int64,Subshell,Subshell,Int64},Vector{Float64}}()
        levels = Level[]
        for  sym  in  relevantSyms
            cache = blockCaches[sym];    idxCsf = cache.idxCsf
            matrix = SelfConsistent.buildCIMatrixEOL(cache, currentOrbitals, grid, nucPot,
                                                      radial1pCache, radial2pCache, vkCache)
            append!( levels, SelfConsistent.diagonalizeBlockEOL(sym, idxCsf, matrix, tempBasis) )
        end
        mp = Basics.sortByEnergy( Multiplet("EOL", levels) )
        return( mp )
    end

    mp           = diagonalizeAllBlocks(orbitals)
    targetLevels = SelfConsistent.selectTargetLevelsEOL(mp, settings.levelSelectionCI)
    previousMc   = [ copy(level.mc)  for level in targetLevels ]
    # THE ROUTE DECIDES whether the off-diagonal CSF-pair terms are kept out of the 1/occ scaling; the global
    # switch of 09-Aug-2026 remains as an override for a caller that names no route.  Leaving it off is what
    # produced the winner-take-all collapse this solver's docstring documents, and it was off by default for a
    # month: measured 09-Sep-2026, turning it on is worth 15.8 mHa on Be 1s^2 2s^2 + 1s^2 2p^2.
    unscaledOff  = settings.scfRoute isa Basics.FockRoute ? settings.scfRoute.unscaledOffDiagonal :
                                                            SelfConsistent.GBL_EOL_UNSCALED_OFFDIAGONAL
    # THE EXPLORATORY ROUTE, and it deliberately uses the GRASP-CONSISTENT SCALING: setcof.f90 divides EVERY
    # contribution -- direct, exchange, one- and two-electron, diagonal and off-diagonal -- by UCF(J), so the
    # off-diagonal terms are NOT left unscaled here.  That is exactly the combination which collapses without
    # stabilization, which is what makes this route the test it is meant to be.
    isStab       = settings.scfRoute isa Basics.StabilizedFockRoute
    if  isStab   unscaledOff = false   end
    oDamp        = Dict{Subshell, Float64}( sh => 0.5  for sh in basis.subshells )
    pedPrev      = Dict{Subshell, Float64}()
    epsPrev      = Dict{Subshell, Float64}()
    if  unscaledOff
        (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels; pairs=:diagonal)
        (coeffs1pOff, coeffs2pOff) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels; pairs=:offdiagonal)
    else
        (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels)
        coeffs2pOff                = Coefficient2p[]
    end
    genOcc = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, targetLevels, basis)

    weights = let  twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
        sumW = sum( twiceJp1(level.J)  for level in targetLevels )
        [ twiceJp1(level.J) / sumW  for level in targetLevels ]
    end
    weightedEnergy = sum( weights[i] * targetLevels[i].energy  for i = 1:length(targetLevels) )

    # (4) Precompute kink-aware Slater-moment tensor caches for every rank that occurs, exactly as
    # solveAverageLevelField does; only the exchange branches of computeTwoElectronV use them
    neededRanks = unique( [ cf.nu for cf in coeffs2p ] )
    if  printout    println(">> [EOL] Precompute kink-aware Slater-moment tensor caches for ranks $(neededRanks) ...")    end
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
        tensorCaches[L] = (cacheLL, cacheLS, cacheSS)
    end

    # THE BEST ITERATE IS KEPT, not merely the last.  This route can walk AWAY from a good orbital set -- that is
    # what the winner-take-all collapse does -- and the CI energy on any orthonormal orbital set is a variational
    # upper bound, so the lowest one seen is the best answer available and is what should be returned.
    bestEnergy   = Inf;    bestOrbitals = deepcopy(orbitals);    bestIteration = 0
    # The b-vectors travel WITH the best orbitals.  They are the same object in two representations, and the
    # orbital energies below are built from the b-vectors, so keeping only one of the pair would pair a
    # reverted orbital set with the b-vectors of a different iteration.
    bestBVectors = deepcopy(bVectors)
    converged    = false;  iterDone     = 0
    for  iter = 1:Basics.maxIterations(settings.scfRoute)
        println("\n> SCF+CI iteration $(iter) [EOL]: ")
        for  level  in  targetLevels
            println("   target level  J=$(level.J)  parity=$(level.parity)  energy=$(level.energy)")
        end

        newBVectors = Dict{Subshell, Vector{Float64}}()
        # Pre-seed with every FROZEN subshell's (fixed) bVector, before any active subshell is refined: this
        # makes Hamiltonian.projectHamiltonian orthogonalize active subshells against frozen same-kappa ones
        # (e.g. a frozen 1s vs. an active, higher-n correlation orbital of the same kappa) and keeps the
        # `count`-based target-eigenvalue-index shift correct, exactly as if the frozen orbitals had already
        # been "processed" this iteration -- which, since they never change, they effectively have.
        processedBVectors = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in settings.frozenSubshells  if  sh in basis.subshells )
        dpm = Dict{Subshell, Float64}()
        # fresh each sweep, for the reason given in solveAverageLevelField: the partner orbitals are fixed
        # within a sweep but not between sweeps, so a cache carried over would serve stale matrices
        directKernels   = Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}}()
        exchangeKernels = Dict{Tuple{Int64,Subshell},Array{Float64,2}}()

        # ORTHY'S ORDERING: spectroscopic orbitals first, correlation orbitals last, so that a weakly occupied
        # orbital is projected against settled ones rather than the other way round.  A stable sort keeps the
        # basis order within each class, so with no correlation orbital the sweep is unchanged.
        sweepSubshells = basis.subshells
        if  isStab
            cut            = settings.scfRoute.correlationCut
            sweepSubshells = sort( collect(basis.subshells), by = sh -> get(genOcc, sh, 0.0) < cut ? 1 : 0 )
        end
        for  subshell  in  sweepSubshells
            if  subshell in settings.frozenSubshells
                # Carry the frozen bVector forward unchanged into newBVectors -- required, since newBVectors
                # is a fresh Dict every iteration and the final newOrbitals tabulation loop below reads
                # bVectors[sh] for EVERY sh in basis.subshells, frozen or not (a KeyError otherwise).
                newBVectors[subshell] = bVectors[subshell]
                continue
            end
            occ = genOcc[subshell]
            # A subshell can legitimately carry ZERO generalized occupation: the EOL functional is built from
            # the target level(s) alone, so a subshell that no target level has any weight on (e.g. 3d when the
            # only target is the J=0+ level of 1s^2 2s^2 + 1s^2 2p^2 + 1s^2 2s3s + 1s^2 2s3d) simply does not
            # enter the energy at all.  The functional is then stationary with respect to that orbital and
            # there is nothing to optimize -- so it is carried forward unchanged, exactly like a frozen one.
            # Refining it regardless used to produce a Fock matrix of NaN (computeFockMatrix divides by
            # occ, and Inf * 0 = NaN), which surfaced far downstream as the thoroughly misleading
            # "Bsplines.findPositiveBranchStart(): no eigenvalue found above the negative-continuum threshold"
            # -- the signature of a missing nuclear well, which was not the problem at all.  AL MEETS THIS TOO
            # and handles it the same way; the claim here that it could not was wrong, because its mean
            # occupation averages over the CSFs of THIS BASIS, and a basis restricted to one symmetry can leave
            # a subshell empty in all of them (Sc+ [Ar] 3d 4s at J = 1 empties 3d_5/2).  Corrected 04-Oct-2026.
            if  abs(occ) < 1.0e-12
                println(">> Subshell $subshell carries zero generalized occupation in the target level(s); " *
                        "the EOL functional does not depend on it, so it is kept unchanged.")
                newBVectors[subshell] = bVectors[subshell]
                continue
            end
            print(">> Refine $subshell orbital with generalized occ = $occ ... ")

            matrix = SelfConsistent.computeFockMatrix(subshell, coeffs2p, bVectors, primitives, nucPot,
                                                              storage, occ, tensorCaches; coeffs2pUnscaled=coeffs2pOff,
                                                              directKernels=directKernels,
                                                              exchangeKernels=exchangeKernels)

            count = Base.count( sh2 -> sh2.kappa == subshell.kappa, keys(processedBVectors) )
            if  count > 0
                matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, processedBVectors)
            end

            wc = Bsplines.diagonalizeLocalMatrix(subshell.kappa, matrix, matrixB, primitives)
            l  = Basics.subshell_l(subshell)
            mm = Bsplines.findPositiveBranchStart(wc.values)
            oldVector = bVectors[subshell]
            # TESTED AND REFUTED (09-Aug-2026): selecting the eigenvector of maximum OVERLAP with the
            # previous orbital instead of the counted index changes nothing here (identical to eight
            # decimals), so the counted index was never the problem -- and it actively harms a correlating
            # orbital, which must be free to change character.  Do not re-propose it.
            ni = mm + subshell.n - l - count - 1
            rawVector = wc.vectors[ni]

            if  transpose(oldVector) * matrixB * rawVector < 0    rawVector = -rawVector    end
            # ADAPTIVE DAMPING, after dampck.f90.  The eigenvalue wc.values[ni] belongs to the orbital just
            # selected and is otherwise discarded; GRASP drives its damping from exactly this quantity.  A SIGN
            # CHANGE between successive relative energy changes is an oscillation and is damped harder, while
            # monotone progress halves the damping so that a settling orbital is taken almost whole.
            damping = 0.5
            if  isStab  &&  settings.scfRoute.adaptiveDamping
                epsNew = wc.values[ni]
                if  haskey(epsPrev, subshell)  &&  abs(epsPrev[subshell]) > 0.
                    ed2 = (epsPrev[subshell] - epsNew) / epsPrev[subshell]
                    if  haskey(pedPrev, subshell)  &&  pedPrev[subshell] * ed2 < -1.0e-4
                        oDamp[subshell] = 0.10 + 0.90 * oDamp[subshell]
                    else
                        oDamp[subshell] = 0.50 * oDamp[subshell]
                    end
                    pedPrev[subshell] = ed2
                end
                epsPrev[subshell] = epsNew
                damping           = oDamp[subshell]
            end
            mixed     = damping * oldVector + (1.0 - damping) * rawVector
            newVector = mixed / sqrt( transpose(mixed) * matrixB * mixed )

            newBVectors[subshell]       = newVector
            processedBVectors[subshell] = newVector

            ovlap = abs( transpose(oldVector) * matrixB * newVector )
            dpm[subshell] = 1.0 - ovlap
            println("     overlap = $ovlap   acc = $(1.0 - ovlap)  ... ")
        end

        bVectors = newBVectors
        newOrbitals = Dict{Subshell, Orbital}()
        # The energy field is 0.0 for the reason given at the sibling solver's final tabulation; the same
        # applies here, and this solver DOES have a one-particle eigenvalue in hand, which is why the choice
        # has to be made deliberately rather than by taking whichever number is nearest.
        for  sh  in  basis.subshells
            newOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives)
        end
        # The EOL driver damps exactly as the AL one does, so it loses same-kappa orthogonality in exactly
        # the same way and needs the same repair.  Measured before this was added: Si^2+ [Ne] 3s^2 + 3p^2
        # gave a worst same-kappa overlap of 2.4e-06 under EOL against 9.3e-10 under AL, and Si^+ reached
        # 6.9e-05.  Wiring the switch into solveAverageLevelField alone was an oversight.
        if  SelfConsistent.GBL_SCF_REORTHONORMALIZE
            (newOrbitals, bVectors) = SelfConsistent.orthonormalizeSameKappa(newOrbitals, bVectors,
                                                        basis.subshells, primitives, matrixB)
        end
        orbitals = newOrbitals

        # Re-diagonalize CI with the refined orbitals; refresh the target level(s)' mixing vectors and energies
        mp           = diagonalizeAllBlocks(orbitals)
        targetLevels = SelfConsistent.selectTargetLevelsEOL(mp, settings.levelSelectionCI)

        # Damping (27-Jul-2026): the outer CI-mixing refresh, feeding straight back into coeffs2p/genOcc every
        # iteration, reproduces the same kind of period-2 oscillation the inner bVector update needed damping
        # for (see solveAverageLevelField/project_df_al_kink_bug.md) -- here in the mixing-vector <->
        # orbital <-> generalized-occupation three-way loop instead of just orbitals <-> orthogonality. Same
        # standard fix: linear mixing of the new and previous mixing vector per target level before use,
        # aligning sign first (a CI eigensolver may return either sign) and renormalizing. The RAW (undamped)
        # targetLevels/energies are still used for reporting and the convergence check -- only the vectors
        # feeding coeffs2p/genOcc are damped.
        dampedLevels = Level[]
        for  (i, level)  in  enumerate(targetLevels)
            newMc = level.mc
            if  transpose(previousMc[i]) * newMc < 0    newMc = -newMc    end
            damping = 0.5
            mixed   = damping * previousMc[i] + (1.0 - damping) * newMc
            mixed   = mixed / sqrt( transpose(mixed) * mixed )
            push!( dampedLevels, Level(level.J, level.M, level.parity, level.index, level.energy,
                                        level.relativeOcc, level.hasStateRep, level.basis, mixed) )
            previousMc[i] = mixed
        end

        if  unscaledOff
            (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, dampedLevels; pairs=:diagonal)
            (coeffs1pOff, coeffs2pOff) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, dampedLevels; pairs=:offdiagonal)
        else
            (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, dampedLevels)
            coeffs2pOff                = Coefficient2p[]
        end
        genOcc = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, dampedLevels, basis)

        weights = let  twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
            sumW = sum( twiceJp1(level.J)  for level in targetLevels )
            [ twiceJp1(level.J) / sumW  for level in targetLevels ]
        end
        newWeightedEnergy = sum( weights[i] * targetLevels[i].energy  for i = 1:length(targetLevels) )

        orbitalConv = maximum( values(dpm) ) < 1.0 ? 1.0 - maximum( values(dpm) ) : 0.0
        energyDiff  = abs( newWeightedEnergy - weightedEnergy )
        println(">> Weighted-average energy = $newWeightedEnergy   orbital-conv = $orbitalConv   " *
                "orbital-acc = $(1.0 - orbitalConv)   energy-diff = $energyDiff")

        weightedEnergy = newWeightedEnergy;    iterDone = iter
        if  newWeightedEnergy < bestEnergy
            bestEnergy = newWeightedEnergy;    bestOrbitals = deepcopy(orbitals);    bestIteration = iter
            bestBVectors = deepcopy(bVectors)
        end
        if  abs(1.0 - orbitalConv) < settings.accuracyScf  &&  energyDiff < settings.accuracyScf
            converged = true;    break
        end
    end

    # THE ROUTE SAYS WHAT IT DID AND WHAT THE ANSWER IS WORTH.  Until 09-Sep-2026 this loop ended in silence,
    # whether it had converged in twelve iterations or run out of budget still moving -- and a quiet wrong answer
    # is the more dangerous of the two failures, because the orbitals it returns look perfectly healthy.
    if  bestEnergy < weightedEnergy - 1.0e-12
        println(">> [EOL-FOCK] the LAST iterate is not the best one: iteration $bestIteration reached " *
                "$bestEnergy against the final $weightedEnergy, so the best is returned.  A Fock iteration that " *
                "moves uphill has been attracted towards the degenerate occ -> 0 fixed point;  treat the result " *
                "with suspicion and compare against Basics.RotationRoute().")
        orbitals = bestOrbitals;    bVectors = bestBVectors
    end
    if  converged
        println(">> [EOL-FOCK] CONVERGED at iteration $iterDone with a weighted-average energy of $bestEnergy.")
    else
        println(">> [EOL-FOCK] STOPPED after $iterDone iterations: the budget of the route was reached with a " *
                "weighted-average energy of $bestEnergy, still moving.  This is NOT convergence;  raise the " *
                "route's maxIterations, or use Basics.RotationRoute(), which descends monotonically.")
        Defaults.warn(AddWarning(), "SelfConsistent.solveOptimizedLevelField(): the EOL field did NOT converge " *
                      "-- $iterDone iterations were exhausted.  The energies are NOT self-consistent.")
    end
    # THE STANDING CAVEAT, printed whether or not the run converged, because it is not a convergence question.
    println(">> [EOL-FOCK] THIS ROUTE IS FAST AND APPROXIMATE.  It converges in tens of iterations where the " *
            "rotation route needs hundreds, and it lands SYSTEMATICALLY ABOVE it: measured 09-Sep-2026 on Be " *
            "1s^2 2s^2 + 1s^2 2p^2, -14.610656 against -14.619515, i.e. 8.86 mHa high.  The cause is known and " *
            "is NOT slow convergence -- the orbital equation, divided by the generalized occupation, has a " *
            "second and degenerate fixed point as that occupation goes to zero.  GRASP avoids it with machinery " *
            "JAC does not have (a node-controlled ODE solution, adaptive damping, Lagrange multipliers, and a " *
            "classification of orbitals into spectroscopic and correlation).  Use Basics.RotationRoute() where " *
            "the number matters.")
    # AND THE SIGNATURE TO CHECK, which is cheap and specific: the two spin-orbit partners of a subshell must
    # agree in mean radius.  On the collapsed answer the Be 2p pair came out at 2.52 and 10.08.
    for  sh  in  basis.subshells
        sh.kappa > 0  ||  continue
        partner = Subshell(sh.n, -sh.kappa - 1)
        haskey(orbitals, partner)  ||  continue
        rA = RadialIntegrals.rkDiagonal(1, orbitals[sh], orbitals[sh], grid)
        rB = RadialIntegrals.rkDiagonal(1, orbitals[partner], orbitals[partner], grid)
        dev = abs(rA - rB) / max(abs(rA), abs(rB))
        # THE THRESHOLD MUST GROW WITH Z, because the partners genuinely split.  A fixed 5 % was used until
        # 09-Sep-2026 and would have fired on EVERY heavy element with perfectly correct orbitals.  Measured on
        # the rotation route, i.e. on orbitals known to be right:  Be at Z = 4 / 10 / 26 / 54 / 92 gives
        # 0.00 / 0.10 / 1.21 / 5.59 / 15.74 %, which for Z >= 26 follows 0.35 (Z alpha)^2.  Carbon does NOT
        # follow it -- 3.06 % at Z = 6 falling to 0.66 % at Z = 26 -- because an open shell splits its partners
        # by their unequal occupations, and that effect is largest where correlation is relatively largest.  So
        # a pure-Z formula cannot serve as the bound; (Z alpha)^2 with a 10 % floor sits about a factor three
        # above every physical value measured, while a genuine collapse (75 % to 154 % in the cases on record)
        # is far above it.
        za    = nuclearModel.Z * Defaults.getDefaults("alpha")
        bound = max(0.10, za*za)
        if  dev > bound
            println(">> [EOL-FOCK] WARNING: the spin-orbit partners $sh and $partner differ in mean radius by " *
                    "$(round(100dev, digits=1)) % ($(round(rA, digits=4)) against $(round(rB, digits=4))), " *
                    "against a bound of $(round(100bound, digits=1)) % for Z = $(nuclearModel.Z).  Some splitting " *
                    "is physical and grows with Z, but not this much;  this is the signature of the collapse " *
                    "described above.")
        end
    end

    # EVERY ORBITAL IS GIVEN ITS DIAGONAL LAGRANGE MULTIPLIER eps_a = <a|F_a|a>, 12-Sep-2026 -- the same
    # quantity, by the same routine, as the rotation solver sets;  see the note there for why it is this and
    # not the one-particle expectation <a|h_D|a>.  This solver has an eigenvalue of F_a in hand during each
    # refinement, and at convergence that eigenvalue and this Rayleigh quotient are the same number;  the
    # quotient is used all the same, so that BOTH solvers report a quantity defined the same way, and so that
    # a run which fell back on an earlier best iterate still reports the multiplier of the orbitals it returns.
    finalTargetLevels = SelfConsistent.selectTargetLevelsEOL(diagonalizeAllBlocks(orbitals), settings.levelSelectionCI)
    if  unscaledOff
        (_, coeffs2pF)    = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, finalTargetLevels; pairs=:diagonal)
        (_, coeffs2pFOff) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, finalTargetLevels; pairs=:offdiagonal)
    else
        (_, coeffs2pF)    = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, finalTargetLevels)
        coeffs2pFOff      = Coefficient2p[]
    end
    genOccF     = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, finalTargetLevels, basis)
    orbEnergies = SelfConsistent.computeOrbitalEnergiesEOL(basis.subshells, bVectors, coeffs2pF, genOccF,
                                                           primitives, grid, nucPot, storage, matrixB;
                                                           coeffs2pUnscaled=coeffs2pFOff)
    for  sh  in  basis.subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, orbEnergies[sh], bVectors[sh], primitives)
    end

    finalBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    multiplet  = Hamiltonian.performCIKinkAware(finalBasis, nuclearModel, grid, settings; printout=printout)
    return( multiplet )
end
