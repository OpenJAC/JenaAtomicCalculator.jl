#
# probeCiCoefficientReuse.jl
#
# A PROBE, not a module:  priority item 17.  `Hamiltonian.setupMatrix` builds a CI matrix by looping over CSF
# pairs and, for every angular coefficient, evaluating
#
#     RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], ...)            (one-particle)
#     InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], ..., cache)              (two-particle)
#
# The two-particle call is memoised -- but its key holds only SUBSHELL LABELS, while the four
# `basis.orbitals[...]` lookups are ARGUMENTS and are therefore evaluated BEFORE the cache is consulted.  With
# a lazy orbital dictionary (PhotoIonizationPseudostates.LazyLadderOrbitals) a lookup can re-materialize a full
# orbital over the whole grid, so on a cache HIT the work is done and thrown away.  The one-particle call has no
# cache at all.
#
# THE QUESTION THIS ANSWERS, and it is combinatorial rather than numerical:  how often does a CI matrix ask for
# the SAME integral again?  That ratio is the fraction of orbital materializations that are pure waste, and it
# needs no orbitals and no radial integrals to measure -- only the angular coefficients.  It is therefore cheap,
# and it can be read on a SMALL grid even though the case that motivated item 17 needs a huge one.
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

function reuse(label, configs, Z)
    grid  = Radial.Grid(true)
    nm    = Nuclear.Model(Z)
    asf   = AsfSettings()
    basis = SelfConsistent.performSCF(configs, nm, grid, asf; printout=false).levels[1].basis
    subshellList = basis.subshells

    syms = Dict{LevelSymmetry,Vector{Int64}}()
    for (i, csf) in enumerate(basis.csfs)
        push!(get!(syms, LevelSymmetry(csf.J, csf.parity), Int64[]), i)
    end

    pairs = 0;  n1 = 0;  n2 = 0
    keys1 = Set{Tuple{Subshell,Subshell}}()
    keys2 = Set{Tuple{Int64,Subshell,Subshell,Subshell,Subshell}}()
    for (sym, idx) in syms
        n = length(idx)
        for r = 1:n,  s = r:n
            pairs += 1
            wa1 = SpinAngular.computeCoefficients(SpinAngular.OneParticleOperator(0, plus),
                                                  basis.csfs[idx[r]], basis.csfs[idx[s]], subshellList)
            wa2 = SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, plus),
                                                  basis.csfs[idx[r]], basis.csfs[idx[s]], subshellList)
            for c in wa1    n1 += 1;  push!(keys1, (c.a, c.b))                      end
            for c in wa2    n2 += 1;  push!(keys2, (c.nu, c.a, c.b, c.c, c.d))      end
        end
    end
    return( (label = label, csfs = length(basis.csfs), subshells = length(subshellList), pairs = pairs,
             n1 = n1, d1 = length(keys1), n2 = n2, d2 = length(keys2)) )
end

cases = [ ("Ne-like Fe, 3 configs",  [Configuration("1s^2 2s^2 2p^6"), Configuration("1s^2 2s^2 2p^5 3s"),
                                      Configuration("1s^2 2s^2 2p^5 3d")], 26.),
          ("Fe 3d^6 4s^2 + 3d^7 4s", [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^6 4s^2"),
                                      Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^7 4s")], 26.),
          ("Si-like Fe, 3 configs",  [Configuration("1s^2 2s^2 2p^6 3s^2 3p^2"),
                                      Configuration("1s^2 2s^2 2p^6 3s 3p^3"),
                                      Configuration("1s^2 2s^2 2p^6 3s^2 3p 3d")], 26.) ]

results = []
for (label, configs, Z) in cases
    @printf("\n>>>> %s\n", label)
    try    push!(results, reuse(label, configs, Z))
    catch err   @printf("     RAISED: %s\n", sprint(showerror, err)[1:min(end,160)])    end
end

println("\n", "="^120)
println("ITEM 17:  how much of setupMatrix's orbital-lookup work is a REPEAT of an integral already computed?")
println("="^120)
@printf("  %-26s %6s %6s %8s  %11s %8s %7s   %11s %8s %7s\n", "case", "CSFs", "subsh", "pairs",
        "1p coeffs", "distinct", "reuse", "2p coeffs", "distinct", "reuse")
println("  ", "-"^118)
for r in results
    @printf("  %-26s %6d %6d %8d  %11d %8d %6.1fx   %11d %8d %6.1fx\n", r.label, r.csfs, r.subshells, r.pairs,
            r.n1, r.d1, r.n1/max(r.d1,1), r.n2, r.d2, r.n2/max(r.d2,1))
end
println("  ", "-"^118)
println("""
  READING.  `reuse` is total coefficients divided by DISTINCT integral keys, i.e. how many times the average
  integral is asked for.  Every one of those requests materializes 2 orbitals (one-particle) or 4 (two-particle)
  BEFORE any cache is consulted, because the lookups are arguments.  With an ordinary Dict of orbitals that is
  cheap; with a lazy, LRU-bounded dictionary over a 850465-point grid it is the dominant cost, and the fraction
  of it that is pure waste is (reuse - 1)/reuse for the two-particle integrals -- and ALL of it for the
  one-particle ones, which have no cache at all.
""")
