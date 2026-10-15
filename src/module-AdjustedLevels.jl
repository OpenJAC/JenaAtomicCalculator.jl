using Printf

"""
`module  JAC.AdjustedLevels`
    ... a submodel of JAC that fits ONE scale factor for the two-body (Coulomb) part of a CI Hamiltonian
        against a handful of measured level energies, re-diagonalizes with the scaled matrix, and returns a
        `FittedMultiplet` -- a distinct type, never a bare `Multiplet`, so that a fitted result can never be
        mistaken for an ab initio one by its type alone.

        WHERE THIS BELONGS IN A CALCULATION, STATED EXPLICITLY SINCE IT IS EASY TO GET WRONG. This operates
        on a `Basis` that already carries CONVERGED ORBITALS from an ordinary SCF (+ RAS, if used) --
        exactly the same `Basis` `Hamiltonian.performCI` itself takes. It does NOT replace or precede the
        SCF step; it replaces the FINAL CI diagonalization at the end of that line, with the SAME orbitals,
        a RESCALED two-body interaction, and a re-fitted set of level energies and mixing coefficients.
        `FittedMultiplet.multiplet` is what should then be handed to whatever property or process
        computation was going to use the ordinary CI result -- `Einstein.computeLines`,
        `IsotopeShift.computeOutcomes`, and so on, exactly as `NistCalibration.CalibratedMultiplet.multiplet`
        already is. No process module depends on this one, nor is this one depended on by any of them --
        the same one-way relationship every module built tonight keeps.

        TWO ROUTES, SELECTED BY WHICH FUNCTION YOU CALL -- never by a flag, the same discipline
        `Hamiltonian.performCI`/`performCIwithFrozenOrbitals`/`performCIKinkAware` already use for choosing
        an algorithm by name. `fitScaleFactor` (this file) fits exactly ONE factor for the WHOLE two-body
        Coulomb matrix element, uniformly across every rank and every pair of orbitals -- the crudest version
        of the idea, chosen because it is the smallest piece that still does the thing that matters (a
        genuine re-diagonalization, so the fit improves the WAVEFUNCTIONS too, not just one level's energy
        the way `NistCalibration` does). `fitOrthogonalParameters` (`module-AdjustedLevels-inc-orthogonal.jl`)
        is the second route: it discovers the SEPARATE integral classes (grouped by multipole rank) actually
        present in a block, fits one parameter per class, and reports which combinations of them the given
        targets can and cannot determine -- the Uylings/Raassen idea, in its most mechanical form. A caller
        who only needs the simple route never has to know the second one exists; the return TYPE
        (`FittedMultiplet` against `FittedMultipletOrthogonal`) says which one produced a given result.

        Also NOT handled: the Breit interaction and QED corrections (pure Coulomb CI only -- the same
        restriction as `settings.eeInteractionCI == CoulombInteraction()`), and fitting ACROSS more than one
        symmetry block at once (the fit works within a single J^P block, against energies measured RELATIVE
        to that block's own lowest level -- there is no cross-block absolute-energy reference here, unlike
        `NistCalibration`'s explicit `groundLevelIndex`). A caller needing either should not read this
        module's numbers as covering that case.

        WHY ONE PARAMETER MAKES THIS TRACTABLE WITHOUT A NEW DEPENDENCY. The expensive part of a CI
        calculation is building the one-body and two-body matrix elements (angular coefficients times radial
        integrals), not diagonalizing the result -- measured elsewhere in this project at roughly 99.6% of
        the cost for a realistic basis. This module builds BOTH matrices ONCE, then searches over the single
        scale factor by just re-combining (oneBody + c*twoBody) and re-diagonalizing for each trial c, which
        is cheap enough that a plain golden-section search (implemented here directly, no external
        optimisation package) over a handful of iterations is entirely sufficient.

        THIS PARALLELS, RATHER THAN CALLS, `Hamiltonian.setupMatrix`. That function sums the one-body and
        two-body contributions into ONE matrix element before returning it, which is exactly the thing this
        module needs to keep SEPARATE in order to scale only the two-body part -- so the same lower-level
        building blocks (`SpinAngular`, `RadialIntegrals.GrantIab`, `InteractionStrength.XL_Coulomb`) are
        reused directly here rather than through `setupMatrix` itself. If that function's own Coulomb-term
        logic changes, this module's `buildOneAndTwoBodyMatrices` may need re-syncing by hand; it is not
        wired to detect that automatically.
"""
module AdjustedLevels

using Printf, LinearAlgebra, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals,
       ..InteractionStrength, ..SpinAngular, ..Hamiltonian


"""
`struct  AdjustedLevels.ScalingParameter`
    ... one fitted scaling parameter. `fitScaleFactor` always returns an Array with exactly ONE of these
        (parameterClass = "Coulomb (uniform)"); `fitOrthogonalParameters` (see this module's
        `-inc-orthogonal` file) returns one PER integral class it discovers, which is why this type -- and
        its sibling `RadialParameter` there -- are named for the ROLE, not the count.

    + parameterClass   ::String    ... which class of integral this factor multiplies.
    + value             ::Float64  ... the fitted factor itself, e.g. 0.82; NOT 1.0 by construction.
"""
struct  ScalingParameter
    parameterClass      ::String
    value               ::Float64
end


"""
`struct  AdjustedLevels.FitTarget`
    ... one measured level this fit was matched against.

    + index            ::Int64    ... position in THIS symmetry block's own energy-sorted spectrum, 1 =
                                        lowest -- NOT a global JAC level index, and NOT comparable across
                                        different symmetry blocks (see this module's own docstring).
    + measuredEnergy   ::Float64  ... the measured excitation energy relative to this SAME block's lowest
                                        level, in cm^-1.
    + fittedEnergy     ::Float64  ... the corresponding excitation energy from the fitted, re-diagonalized
                                        spectrum, in cm^-1, for direct comparison.
"""
struct  FitTarget
    index               ::Int64
    measuredEnergy       ::Float64
    fittedEnergy         ::Float64
end


"""
`struct  AdjustedLevels.FittedMultiplet`
    ... the result of `AdjustedLevels.fitScaleFactor`. See this module's own docstring for what this is, and
        -- at least as importantly -- is NOT.

    + multiplet       ::ManyElectron.Multiplet     ... the single-symmetry-block multiplet from the FITTED,
                                                         re-diagonalized Hamiltonian; same CSFs/orbitals as
                                                         the input Basis.
    + scaleFactors    ::Array{ScalingParameter,1}   ... the fitted scaling parameter(s) -- never empty.
    + fitTargets       ::Array{FitTarget,1}          ... the measured levels this fit was matched against.
    + rmsDeviation      ::Float64                    ... root-mean-square deviation, fitted vs measured, over
                                                           the fit targets, in cm^-1.
    + source            ::String                    ... a dated citation for the measured data used.
"""
struct  FittedMultiplet
    multiplet            ::ManyElectron.Multiplet
    scaleFactors         ::Array{ScalingParameter,1}
    fitTargets           ::Array{FitTarget,1}
    rmsDeviation         ::Float64
    source               ::String
end


"""
`Base.show(io::IO, fm::AdjustedLevels.FittedMultiplet)`
    ... ALWAYS prints the "fitted, not ab initio" disclosure first, even on an incidental `println(fm)` or
        `@show fm` -- the one design choice this module makes beyond what `NistCalibration` does, because a
        skippable disclosure is not loud enough for a result built from adjustable parameters. Nothing is
        returned.
"""
function  Base.show(io::IO, fm::FittedMultiplet)
    println(io, "FittedMultiplet -- FITTED, NOT AB INITIO.")
    println(io, "  source: $(fm.source)")
    for  sf in fm.scaleFactors
        println(io, "  scale factor [$(sf.parameterClass)] = $(sf.value)")
    end
    println(io, "  fit against $(length(fm.fitTargets)) measured level(s), rms deviation = $(fm.rmsDeviation) cm^-1")
    println(io, "  multiplet: $(fm.multiplet.name), $(length(fm.multiplet.levels)) level(s)")
end


"""
`AdjustedLevels.buildOneAndTwoBodyMatrices(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid)`
    ... to build the one-body and two-body (PURE COULOMB) CI matrices for the given symmetry block,
        separately -- see this module's own docstring for why this parallels rather than calls
        `Hamiltonian.setupMatrix`. A (idxCsf, oneBody, twoBody) tuple is returned, with idxCsf::Array{Int64,1}
        the basis.csfs indices of this block and oneBody, twoBody::Array{Float64,2} upper-triangular (the
        same convention `Hamiltonian.setupMatrix` uses, relying on `Basics.diagonalize`'s `Symmetric` wrapper
        to read only that triangle).
"""
function  buildOneAndTwoBodyMatrices(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid)
    idxCsf = [idx for idx = 1:length(basis.csfs) if basis.csfs[idx].J == sym.J && basis.csfs[idx].parity == sym.parity]
    n = length(idxCsf)
    potential = Nuclear.nuclearPotential(nm, grid)
    oneBody = zeros(Float64, n, n);   twoBody = zeros(Float64, n, n)

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
            me2 = 0.0
            for  coeff in waG2
                me2 = me2 + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                               basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
            end
            oneBody[r,s] = me1;   twoBody[r,s] = me2
        end
    end

    return( (idxCsf, oneBody, twoBody) )
end


"""
`AdjustedLevels.goldenSectionMinimize(f::Function, a::Float64, b::Float64; tolerance::Float64=1.0e-4, maxIterations::Int64=60)`
    ... a plain, dependency-free golden-section search for the minimum of `f` on `[a,b]`, assumed unimodal
        there -- sufficient for the smooth, single-parameter scan this module needs, see its own docstring.
        A (xBest, fBest)::Tuple{Float64,Float64} is returned.
"""
function  goldenSectionMinimize(f::Function, a::Float64, b::Float64; tolerance::Float64=1.0e-4, maxIterations::Int64=60)
    gr = (sqrt(5.0) - 1.0) / 2.0
    c = b - gr*(b-a);   d = a + gr*(b-a)
    fc = f(c);   fd = f(d)
    for  _ = 1:maxIterations
        (b - a) < tolerance   &&   break
        if  fc < fd
            b = d;   d = c;   fd = fc;   c = b - gr*(b-a);   fc = f(c)
        else
            a = c;   c = d;   fc = fd;   d = a + gr*(b-a);   fd = f(d)
        end
    end
    xBest = fc < fd ? c : d

    return( (xBest, f(xBest)) )
end


"""
`AdjustedLevels.fitScaleFactor(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                targets::Array{Tuple{Int64,Float64},1}, source::String;
                                scaleBounds::Tuple{Float64,Float64}=(0.5,1.2))`
    ... to fit ONE Coulomb scale factor within `scaleBounds` so that the re-diagonalized spectrum of `sym`'s
        block matches `targets` (pairs of (position in the sorted spectrum, measured energy relative to the
        block's own lowest level, in cm^-1)) as closely as possible in the least-squares sense. A
        fm::AdjustedLevels.FittedMultiplet is returned.
"""
function  fitScaleFactor(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                          targets::Array{Tuple{Int64,Float64},1}, source::String;
                          scaleBounds::Tuple{Float64,Float64}=(0.5,1.2))
    isempty(targets)   &&   error("AdjustedLevels.fitScaleFactor(): need at least one target level to fit against.")
    idxCsf, oneBody, twoBody = buildOneAndTwoBodyMatrices(sym, basis, nm, grid)
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)

    function  relativeEnergiesKayser(c::Float64)
        eigen = Hamiltonian.diagonalizeCiMatrix(oneBody .+ c .* twoBody, LevelSelection())
        energies = sort(eigen.values)
        return( (energies .- energies[1]) ./ auPerKayser )
    end

    function  rms(c::Float64)
        rel = relativeEnergiesKayser(c)
        s = 0.0
        for  (idx, measured) in targets   s = s + (rel[idx] - measured)^2   end

        return( sqrt(s/length(targets)) )
    end

    cBest, rmsBest = goldenSectionMinimize(rms, scaleBounds[1], scaleBounds[2])

    eigen  = Hamiltonian.diagonalizeCiMatrix(oneBody .+ cBest .* twoBody, LevelSelection())
    order  = sortperm(eigen.values)
    energiesAU  = eigen.values[order]
    relKayser   = (energiesAU .- energiesAU[1]) ./ auPerKayser

    levels = Level[]
    for  (pos, ev) in enumerate(order)
        vector = zeros( length(basis.csfs) );   ns = 0
        for  r = 1:length(basis.csfs)
            if  LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = eigen.vectors[ev][ns]   end
        end
        push!( levels, Level(sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, pos, energiesAU[pos], 0., true, basis, vector) )
    end

    fitTargets = [ FitTarget(idx, measured, relKayser[idx]) for (idx, measured) in targets ]

    return( FittedMultiplet(Multiplet(string(sym)*"+ (fitted)", levels), [ScalingParameter("Coulomb (uniform)", cBest)],
                             fitTargets, rmsBest, source) )
end


"""
`AdjustedLevels.mergeIntoMultiplet(fm::AdjustedLevels.FittedMultiplet, abInitioMultiplet::ManyElectron.Multiplet)`
    ... to combine a FittedMultiplet's single, re-fitted symmetry block with the REST of an ordinary, ab initio
        multiplet -- typically the very same `Hamiltonian.performCI` result `fitScaleFactor`'s own `basis` came
        from, which carries every OTHER symmetry block a fit never touches, such as a ground state of different
        J^P. This is the step named in this module's own docstring as happening at the end of an SCF+RAS line:
        the levels of `abInitioMultiplet` belonging to `fm`'s own symmetry are replaced outright by
        `fm.multiplet.levels`, every other symmetry block is kept unchanged, and the result is re-sorted by
        energy exactly as `Hamiltonian.performCI` itself does after merging its own per-block multiplets -- so
        that the combined multiplet can be handed to `Einstein.computeLines` or any other property/process
        computation needing levels from more than one symmetry. A mp::ManyElectron.Multiplet is returned.
"""
function  mergeIntoMultiplet(fm::FittedMultiplet, abInitioMultiplet::ManyElectron.Multiplet)
    isempty(fm.multiplet.levels)   &&   error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultiplet has no levels.")
    sym = LevelSymmetry(fm.multiplet.levels[1].J, fm.multiplet.levels[1].parity)
    otherLevels = [lev for lev in abInitioMultiplet.levels if LevelSymmetry(lev.J, lev.parity) != sym]

    length(otherLevels) + length(fm.multiplet.levels) != length(abInitioMultiplet.levels)   &&
        error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultiplet's block has a different number of levels than " *
              "the matching block in abInitioMultiplet -- they must come from the same Basis.")

    allLevels = vcat(otherLevels, fm.multiplet.levels)
    mp = Multiplet(abInitioMultiplet.name * "+fitted", allLevels)

    return( Basics.sortByEnergy(mp) )
end


include("module-AdjustedLevels-inc-orthogonal.jl")

end # module
