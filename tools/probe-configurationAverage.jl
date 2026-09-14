#
# probe-configurationAverage.jl   --   the CLOSED FORM for the average energy of a configuration, checked against
#                                      the exact (2J+1)-weighted trace.
#
# PRIORITY ITEM 22.  `Basics.determineMeanEnergy` has two routes today: one diagonalises the configuration's
# multiplet, and one (14-Sep-2026, `4bb483c`) takes the (2J+1)-weighted trace of the CI matrix, which is EXACT --
# a similarity transformation preserves the trace and the matrix is block-diagonal in J^P -- and removes the
# eigenvalue work.  Neither removes the CSF ENUMERATION.  Slater's average energy of a configuration is O(1) in
# the number of subshell PAIRS and touches no CSF at all, and that is what this file works out.
#
# WHY IT MUST BE DERIVED RATHER THAN TRANSCRIBED.  JAC's orbitals are jj-coupled while a `Configuration` is not,
# so the textbook non-relativistic E_av does not carry over.  The relativistic pieces used here are
#
#     one-electron          SUM_a  w_a I(a)
#     same subshell a       w_a(w_a-1)/2 * [ F^0(aa) - SUM_{k=2,4,..,2j_a} <a||C^k||a>^2 / (2j_a (2j_a+1)) F^k(aa) ]
#     different a, b        w_a w_b      * [ F^0(ab) -       SUM_k <a||C^k||b>^2 / ((2j_a+1)(2j_b+1)) G^k(ab) ]
#
# BOTH EXCHANGE PREFACTORS ARE UNITY, and that was established rather than assumed.  A first attempt carried the
# textbook non-relativistic 1/2 on the cross term and failed: measured 14-Sep-2026 on closed shells, where the
# trace is exact, 1s^2 2s^2 (which has NO same-subshell term, so it isolates the cross one) gave residual
# -0.05062436 against S_cross 0.05062436, i.e. a factor of exactly 1.0000;  with that fixed, 1s^2 2s^2 2p^6 and
# [Ar] give the same-subshell factor as 1.0000 and 0.99999883.  So there was ONE error, the spurious 1/2, and
# nothing was tuned -- two constants that both come out unity to seven figures are not a fit.
#
# THE INTEGRAL MUST BE THE KINK-AWARE ONE.  With the plain `SlaterRk`, 1s^2 -- which has no k > 0 term at all and
# must therefore be exact -- came out 1.1e-04 off, a relative 4.9e-05.  That is quadrature, not physics: the CI
# path the trace uses builds its Slater integrals with `SlaterRkKinkAware`, and comparing against it with a
# different quadrature measures the quadrature.  With it, 1s^2 agrees to 1.8e-15.
#
# with F^k(ab) = R^k(abab) and G^k(ab) = R^k(abba).  A j = 1/2 subshell admits no even k > 0, so 1s^2 reduces to
# I(1s) x 2 + F^0(1s,1s), which is the right answer and the cheapest check that the prefactors are not nonsense.
#
# AND THE NON-RELATIVISTIC CONFIGURATION IS A WEIGHTED SUM OVER THE RELATIVISTIC ONES.  `3d^2` is not one jj
# configuration but every way of splitting two electrons over 3d_3/2 and 3d_5/2, and the average over the
# non-relativistic configuration weights each by its number of magnetic substates, PROD_a binomial(2j_a+1, w_a).
# That weight is exactly the denominator the trace identity uses -- SUM_csfs (2J_r+1) over one relativistic
# configuration IS its determinant count -- so the two routes are being compared on the same average and any
# disagreement is the formula's, not the weighting's.
#
# THE TEST IS COLUMN 3.  The trace route is exact, so the closed form must reproduce it to round-off on a closed
# shell and to the trace identity on an open one.  A systematic ratio rather than noise means a prefactor is wrong.
#
using JenaAtomicCalculator, Printf
const B  = JenaAtomicCalculator.Basics
const RI = JenaAtomicCalculator.RadialIntegrals
const AM = JenaAtomicCalculator.AngularMomentum

"binomial weight of one relativistic configuration = its number of magnetic substates"
function detCount(rconf::ConfigurationR)
    w = 1.0
    for (sh, occ) in rconf.subshells
        occ == 0  &&  continue
        w *= binomial(Basics.subshell_2j(sh) + 1, occ)
    end
    return( w )
end

# THE RADIAL INTEGRALS ARE CACHED ACROSS THE RELATIVISTIC CONFIGURATIONS, and without that the closed form is
# SLOWER than the trace it is meant to replace -- measured 14-Sep-2026: 3d^5 took 0.357 s against the trace's
# 0.12 s, because `3d^5` splits many ways over 3d_3/2 and 3d_5/2 and every one of them rebuilt the same F^k and
# G^k.  The integrals depend only on (k, subshell, subshell), never on the occupations, so one Dict over the
# whole configuration removes the repetition entirely.
const RKCACHE = Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}()
function rk(k, a::Subshell, b::Subshell, c::Subshell, d::Subshell, orbitals, grid)
    key = (k, a, b, c, d)
    haskey(RKCACHE, key)  &&  return( RKCACHE[key] )
    v = RI.SlaterRkKinkAware(k, orbitals[a], orbitals[b], orbitals[c], orbitals[d], grid)
    RKCACHE[key] = v
    return( v )
end
const I1CACHE = Dict{Subshell,Float64}()
function i1(a::Subshell, orbitals, grid, nucPot)
    haskey(I1CACHE, a)  &&  return( I1CACHE[a] )
    v = RI.GrantIab(orbitals[a], orbitals[a], grid, nucPot);    I1CACHE[a] = v
    return( v )
end


"the closed-form average energy of ONE relativistic configuration"
function eavRelativistic(rconf::ConfigurationR, orbitals, grid, nucPot)
    shs = [ sh for (sh, occ) in rconf.subshells  if occ > 0 ]
    sort!(shs, by = s -> (s.n, s.kappa))
    E = 0.0
    for sh in shs
        w = rconf.subshells[sh]
        E += w * i1(sh, orbitals, grid, nucPot)
    end
    for (ia, sha) in enumerate(shs)
        wa = rconf.subshells[sha];    ja2 = Basics.subshell_2j(sha);    oa = orbitals[sha]
        if  wa >= 2
            F0 = rk(0, sha, sha, sha, sha, orbitals, grid);    s = 0.0
            for k = 2:2:ja2
                c = AM.CL_reduced_me(sha, k, sha)
                c == 0.  &&  continue
                s += c^2 / (ja2 * (ja2 + 1)) * rk(k, sha, sha, sha, sha, orbitals, grid)
            end
            E += wa*(wa-1)/2 * (F0 - s)
        end
        for  ib = (ia+1):length(shs)
            shb = shs[ib];   wb = rconf.subshells[shb];   jb2 = Basics.subshell_2j(shb);   ob = orbitals[shb]
            F0 = rk(0, sha, shb, sha, shb, orbitals, grid);    s = 0.0
            for k = 0:((ja2 + jb2) ÷ 2)
                c = AM.CL_reduced_me(sha, k, shb)
                c == 0.  &&  continue
                s += c^2 / ((ja2 + 1) * (jb2 + 1)) * rk(k, sha, shb, shb, sha, orbitals, grid)
            end
            E += wa * wb * (F0 - s)
        end
    end
    return( E )
end

"the closed-form average over a NON-relativistic configuration"
function eavClosedForm(conf::Configuration, orbitals, nm, grid)
    nucPot  = Nuclear.nuclearPotential(nm, grid)
    empty!(RKCACHE);    empty!(I1CACHE)          ## the orbitals change between cases, so the cache is per call
    rconfs  = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
    num = 0.0;   den = 0.0
    for rc in rconfs
        w = detCount(rc);    num += w * eavRelativistic(rc, orbitals, grid, nucPot);    den += w
    end
    return( num / den )
end

# THE LAST FOUR ARE THE POINT.  The closed form never touches a CSF, so its cost is set by the number of subshell
# PAIRS and is flat in the CSF count, while the trace must enumerate every CSF.  On the small cases that shows as
# little -- 3d^5 has 37 CSFs -- so the open f shells, where one configuration carries hundreds, are where the
# claim is actually tested.  The CSF count is printed beside the timings so the two can be read together.
cases = [("closed  1s^2",                  4.0,  "1s^2"),
         ("closed  1s^2 2s^2 2p^6",       10.0,  "1s^2 2s^2 2p^6"),
         ("open    1s^2 2s 2p",            6.0,  "1s^2 2s 2p"),
         ("open    1s^2 2s^2 2p^3",        7.0,  "1s^2 2s^2 2p^3"),
         ("open    [Ar] 3d^2",            22.0,  "1s^2 2s^2 2p^6 3s^2 3p^6 3d^2"),
         ("open    [Ar] 3d^5",            25.0,  "1s^2 2s^2 2p^6 3s^2 3p^6 3d^5"),
         ("open    [Ar] 3d^7",            27.0,  "1s^2 2s^2 2p^6 3s^2 3p^6 3d^7"),
         ("open    [Xe] 4f^5",            63.0,  "1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^5"),
         ("open    [Xe] 4f^7",            65.0,  "1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^7")]

@printf("\n%-26s %7s %18s %18s %12s %9s %9s %7s\n",
        "configuration", "CSFs", "trace (exact)", "closed form", "difference", "trace s", "closed s", "gain")
println("-"^116)
for (tag, Z, cstr) in cases
    conf = Configuration(cstr);   nm = Nuclear.Model(Z)
    grid = Basics.recommendedGrid([conf], nm; printout=false)
    set  = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction(),
                                      eeInteractionCI = CoulombInteraction(), gridStopper = false)
    mp   = redirect_stdout(devnull) do
               SelfConsistent.performSCF([conf], nm, grid, set; printout=false)
           end
    orbs = mp.levels[1].basis.orbitals
    t1 = @elapsed (eTrace = redirect_stdout(devnull) do
                                Basics.determineMeanEnergy(conf, orbs, nm, grid, set)
                            end)
    t2 = @elapsed (eClosed = eavClosedForm(conf, orbs, nm, grid))
    nCsf = length(mp.levels[1].basis.csfs)
    @printf("%-26s %7d %18.9f %18.9f %12.2e %9.2f %9.3f %6.1fx\n",
            tag, nCsf, eTrace, eClosed, abs(eTrace-eClosed), t1, t2, t1/max(t2,1e-9))
    flush(stdout)
end
println("\nThe trace route is EXACT, so column 3 is the verdict:  round-off means the closed form is right;  a")
println("systematic size, and above all one that tracks the open-shell occupation, means a prefactor is wrong.")
