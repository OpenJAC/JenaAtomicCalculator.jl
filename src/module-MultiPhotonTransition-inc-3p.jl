# Three-photon excitation and decay between bound levels.
# WHAT IS AND IS NOT IMPLEMENTED HERE, 08-Aug-2026, because the two halves differ:
#   * THREE-PHOTON ABSORPTION is implemented, in an elementary formulation -- see the second header, further
#     down, for exactly what "elementary" buys and what it gives up.
#   * THREE-PHOTON EMISSION is NOT: its computeLines still refuses, with a message naming what would have to be
#     built. What exists for it is the piece that could be finished and VERIFIED on its own, the ENERGY SHARINGS
#     ON THE TWO-DIMENSIONAL SIMPLEX.
# The asymmetry is not an oversight. Absorption fixes the three photon energies, so it needs no sharings at all;
# emission fixes only their SUM, so it needs the simplex, and it is that quadrature -- not the amplitude -- that
# was the self-contained piece worth doing first.
# WHY THAT PIECE FIRST. For two photons the sharings are a Gauss-Legendre line: omega1 in [0, E], omega2 fixed by
# energy conservation. For three photons the domain is the triangle
#     omega1 + omega2 + omega3 = E,      omega1, omega2, omega3 > 0
# and the rate is differential in TWO of them. Getting that quadrature right is a self-contained problem with an
# EXACT answer to check against, which is what makes it worth doing before the third-order amplitude: the moments
# of a simplex are known in closed form,
#     Int omega1^a omega2^b omega3^c domega1 domega2  =  E^(a+b+c+2) * a! b! c! / (a+b+c+2)!
# so a = b = c = 0 must give E^2/2, (1,0,0) must give E^3/6, and (1,1,1) must give E^5/120. Three sum rules of
# increasing polynomial degree, all parameter-free. If the weights are wrong, the eventual three-photon total
# rate is wrong by a factor that no amount of checking the amplitude would ever reveal -- exactly the class of
# silent normalisation error that has cost this module time before.
# THE TRANSFORM is the standard collapsed-coordinate (Duffy) map of the square onto the triangle,
#     omega1 = E u,     omega2 = E (1-u) v,     omega3 = E (1-u)(1-v),     u, v in [0,1]
# with Jacobian domega1 domega2 = E^2 (1-u) du dv. A tensor Gauss-Legendre rule in (u,v) then integrates any
# polynomial exactly up to the order of the rule, including the (1-u) weight. `noSharings` is the number of
# points PER DIRECTION, so the grid carries noSharings^2 sharings -- stated because it is the one place a user
# could reasonably expect it to mean the total.
# ONE HONEST LIMITATION, stated rather than discovered later: the Duffy grid is EXACT but not PERMUTATION
# SYMMETRIC -- omega1 plays a distinguished role in the transform. For the integrated rate that is irrelevant
# (exactness is exactness), but a three-photon spectrum plotted on these nodes will not LOOK symmetric under
# exchanging the photons, even though the underlying function is. The symmetry check that served the two-photon
# case so well (blocker A1 was found through it) therefore has to be applied to the function, evaluated at
# permuted points, and not to the tabulated node values.


"""
`struct  MultiPhotonTransition.Sharing_3p`
    ... defines a type for one energy sharing of the total transition energy among THREE photons.

    + omega1           ::Float64     ... Energy of photon 1.
    + omega2           ::Float64     ... Energy of photon 2.
    + omega3           ::Float64     ... Energy of photon 3.
    + weight           ::Float64     ... Quadrature weight of this sharing for energy-integrated quantities,
                                         i.e. sum_i weight_i f_i approximates the integral of f over the simplex
                                         with respect to domega1 domega2.
"""
struct  Sharing_3p
    omega1             ::Float64
    omega2             ::Float64
    omega3             ::Float64
    weight             ::Float64
end


# `Base.show(io::IO, sharing::MultiPhotonTransition.Sharing_3p)`
#   ... prepares a proper printout of the variable sharing::MultiPhotonTransition.Sharing_3p.
function Base.show(io::IO, sharing::MultiPhotonTransition.Sharing_3p)
    println(io, "omega1:            $(sharing.omega1)  ")
    println(io, "omega2:            $(sharing.omega2)  ")
    println(io, "omega3:            $(sharing.omega3)  ")
    println(io, "weight:            $(sharing.weight)  ")
end


"""
`MultiPhotonTransition.determineSharings_3p(energy::Float64, noSharings::Int64)`
    ... to determine the energy sharings of the total energy among three photons, placed on a tensor
        Gauss-Legendre grid over the two-dimensional simplex; an Array{MultiPhotonTransition.Sharing_3p,1} with
        noSharings^2 entries is returned.

        The weights integrate with respect to domega1 domega2 and already contain the Jacobian of the
        collapsed-coordinate map; see the header of this file and `checkSharings_3p`, which verifies them
        against the exact simplex moments rather than asserting them.
"""
function determineSharings_3p(energy::Float64, noSharings::Int64)
    if  noSharings < 1    error("MultiPhotonTransition: noSharings must be >= 1, not $noSharings.")    end
    if  energy <= 0.      error("MultiPhotonTransition: the total energy must be positive, not $energy.")    end
    sharings = MultiPhotonTransition.Sharing_3p[]
    xx, ww   = GaussQuadrature.legendre(noSharings)
    for  (iu, xu) in enumerate(xx)
        u  = (xu + 1.) / 2;    wu = ww[iu] / 2       ## map [-1,1] -> [0,1]
        for  (iv, xv) in enumerate(xx)
            v  = (xv + 1.) / 2;    wv = ww[iv] / 2
            omega1 = energy * u
            omega2 = energy * (1. - u) * v
            omega3 = energy * (1. - u) * (1. - v)
            # the Jacobian domega1 domega2 = E^2 (1-u) du dv
            weight = wu * wv * energy^2 * (1. - u)
            push!( sharings, MultiPhotonTransition.Sharing_3p(omega1, omega2, omega3, weight) )
        end
    end
    return( sharings )
end


"""
`MultiPhotonTransition.checkSharings_3p(energy::Float64, noSharings::Int64)`
    ... to verify the simplex quadrature against the exact moments of the simplex; an
        errors::Array{Tuple{String,Float64,Float64,Float64},1} of (moment, exact, quadrature, relative error) is
        returned.

        THIS IS THE SUM RULE, and it is what makes the sharings a finished piece of work rather than plausible
        code. The moments of the two-dimensional simplex are known exactly,

            Int omega1^a omega2^b omega3^c domega1 domega2  =  E^(a+b+c+2) a! b! c! / (a+b+c+2)!

        so (0,0,0) -> E^2/2, (1,0,0) -> E^3/6, (1,1,1) -> E^5/120. The three cases probe polynomial degrees 0, 1
        and 3, and the last one is the sharp one: it is degree 3 in the collapsed coordinates and would expose a
        missing or misplaced Jacobian, which the constant test alone cannot.
"""
function checkSharings_3p(energy::Float64, noSharings::Int64)
    sharings = MultiPhotonTransition.determineSharings_3p(energy, noSharings)
    results  = Tuple{String,Float64,Float64,Float64}[]
    for  (name, a, b, c, exact) in [("1",             0, 0, 0, energy^2 / 2),
                                    ("omega1",        1, 0, 0, energy^3 / 6),
                                    ("omega1 omega2", 1, 1, 0, energy^4 / 24),
                                    ("w1 w2 w3",      1, 1, 1, energy^5 / 120)]
        quad = 0.
        for  sh in sharings
            quad = quad + sh.weight * sh.omega1^a * sh.omega2^b * sh.omega3^c
        end
        push!( results, (name, exact, quad, abs(quad - exact) / abs(exact)) )
    end
    return( results )
end


"""
`MultiPhotonTransition.displaySharings_3p(stream::IO, energy::Float64, noSharings::Int64)`
    ... to display the three-photon energy sharings together with the sum-rule check of their weights; nothing
        is returned.
"""
function  displaySharings_3p(stream::IO, energy::Float64, noSharings::Int64)
    sharings = MultiPhotonTransition.determineSharings_3p(energy, noSharings)
    nx = 76
    println(stream, " ")
    println(stream, "  Three-photon energy sharings on the simplex omega1 + omega2 + omega3 = " *
            "$(round(Defaults.convertUnits("energy: from atomic", energy), sigdigits=6)) " *
            TableStrings.inUnits("energy") * "  ($(noSharings) points per direction, " *
            "$(length(sharings)) in total):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "        omega1          omega2          omega3            weight")
    println(stream, "        " * TableStrings.inUnits("energy") * "            " *
                    TableStrings.inUnits("energy") * "            " * TableStrings.inUnits("energy"))
    println(stream, "  ", TableStrings.hLine(nx))
    for  sh in sharings
        println(stream, "     " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", sh.omega1)) *
                "    " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", sh.omega2)) *
                "    " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", sh.omega3)) *
                "    " * @sprintf("%.6e", sh.weight))
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, " ")
    println(stream, "  Sum rule -- the quadrature against the EXACT moments of the simplex:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "     moment              exact           quadrature      relative error")
    println(stream, "  ", TableStrings.hLine(nx))
    for  (name, exact, quad, err) in MultiPhotonTransition.checkSharings_3p(energy, noSharings)
        println(stream, "     " * TableStrings.flushleft(16, name; na=2) * @sprintf("%.8e", exact) * "    " *
                @sprintf("%.8e", quad) * "    " * @sprintf("%.3e", err))
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">>> All four moments must be reproduced to machine precision; the (1,1,1) moment is the " *
            "one that\n>>> would expose a missing Jacobian, which the constant moment alone cannot.")
    return( nothing )
end


"""
`MultiPhotonTransition.computeLines(scheme::ThreePhotonEmissionScheme, finalMultiplet::Multiplet,
                            initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings; output=true)`
    ... to compute three-photon emission lines. The AMPLITUDE is not implemented and an informative error is
        raised -- except under `calcOverview`, which displays the piece that IS implemented and verified, namely
        the energy sharings on the simplex and their sum-rule check.
"""
function  computeLines(scheme::ThreePhotonEmissionScheme, finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                       grid::Radial.Grid, settings::MultiPhotonTransition.Settings; output=true)
    # calcOverview shows what exists rather than refusing outright. The sharings need nothing but the transition
    # energy, so they can be inspected and checked long before any amplitude is written -- which is the point of
    # having built them first.
    if  settings.calcOverview
        println("")
        printstyled("MultiPhotonTransition.computeLines(::ThreePhotonEmissionScheme): only the energy sharings " *
                    "are implemented; showing them. \n", color=:light_green)
        println("")
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    energy = iLevel.energy - fLevel.energy + settings.photonEnergyShift
                    if  energy <= 0.    continue    end
                    println("  Line $(iLevel.index) - $(fLevel.index):")
                    MultiPhotonTransition.displaySharings_3p(stdout, energy, scheme.noSharings)
                end
            end
        end
        println("\n>>> Overview only; the third-order amplitude is not implemented, so no rate was computed.")
        if  output    return( MultiPhotonTransition.Sharing_3p[] )
        else          return( nothing )
        end
    end
    error("\n\nMultiPhotonTransition.computeLines():  STOP -- ThreePhotonEmissionScheme is not yet implemented.\n" *
          ">>> What is missing, in the order it would have to be built:\n"                                        *
          "    (1) a THIRD-order amplitude, i.e. a double sum over two sets of intermediate states\n"             *
          "        <f|O(mp3)|nu2> <nu2|O(mp2)|nu1> <nu1|O(mp1)|i> / [(E_i - w1 - E_nu1)(E_i - w1 - w2 - E_nu2)]\n"*
          "        with all 3! = 6 orderings of the photons, against 2! = 2 for the two-photon case;\n"           *
          "    (3) the coupling of THREE photon multipoles, so the rank K is no longer fixed by\n"                *
          "        oplus(J_f, J_i) alone and an intermediate coupling rank must be summed over.\n"                *
          ">>> WHAT IS ALREADY DONE (08-Aug-2026): step (2), the energy sharings on the two-dimensional simplex,\n"*
          "    verified against the exact simplex moments. Run the same scheme with calcOverview = true to see\n" *
          "    them and their sum-rule check, or call MultiPhotonTransition.determineSharings_3p(energy, n).\n"   *
          ">>> Use TwoPhotonEmissionScheme() meanwhile.\n")
end




# ---------------------------------------------------------------------------------------------------------------
# THREE-PHOTON ABSORPTION, in an ELEMENTARY formulation.       Added 08-Aug-2026.
# WHAT "ELEMENTARY" MEANS HERE, stated first so that nothing below is mistaken for more than it is:
#   * THREE BEAMS, ALL LINEARLY POLARIZED ALONG THE SAME AXIS. One geometry, no polarization decomposition, no
#     Stokes parameters. With every polarization along z only the q = 0 spherical component contributes, the
#     magnetic quantum number is conserved through the whole chain, and the four-fold m-sum collapses to a
#     single M. That collapse is what makes this case worth doing first.
#   * A STRENGTH, NOT A CROSS SECTION.  S^(3) = 1/(2J_i+1) sum_M |A(M)|^2 in atomic units. A generalized
#     three-photon cross section needs a normalisation of order F^3; the module's TWO-photon absorption
#     normalisation has never been derived either, so inventing a three-photon one would add a second
#     undetermined constant wearing the units of a measured quantity.
#   * NO POLARIZATION OBSERVABLES and no rank-K decomposition.
# WHY THE m-SUM AND NOT A COUPLED-TENSOR FORM. The natural-looking route is to couple the three photon
# multipoles to a total rank K through an intermediate rank k, as the two-photon code couples two. It is also
# the route on which this implementation would most plausibly have been WRONG: the six time orderings then have
# to be re-expressed in one common coupling order, which brings in recoupling coefficients (6-j, and 9-j for
# unequal multipoles) whose phases are easy to get confidently wrong and hard to test. Summing over magnetic
# quantum numbers with one Wigner-Eckart 3-j per step needs no recoupling at all: each ordering is a plain
# product of three 3-j symbols and three reduced matrix elements. It costs an m-sum -- which for parallel
# linear polarization is a single loop -- and buys transparency.
# NOTE WHAT IS THEREFORE *NOT* A TEST HERE. For two photons, invariance of the result under exchanging the two
# colours was a genuine check, because the two orderings were combined with a relative phase that could be (and
# was) wrong. Here all six orderings are summed explicitly with no relative phase to get wrong, so permuting the
# three colours is invariant BY CONSTRUCTION and proves nothing. The checks that do have teeth are the parity
# selection rule (an exact zero), the count of six orderings in the monochromatic limit, and the hydrogenic
# Z-scaling; see examples/example-Du.jl.
# ---------------------------------------------------------------------------------------------------------------


"""
`struct  MultiPhotonTransition.Line_3pAbsorption`
    ... defines a type for a three-photon absorption line driven by three parallel linearly-polarized beams.

    + initialLevel     ::Level        ... initial-(state) level
    + finalLevel       ::Level        ... final-(state) level
    + omega1           ::Float64      ... energy of the first photon.
    + omega2           ::Float64      ... energy of the second photon.
    + omega3           ::Float64      ... energy of the third photon.
    + strength         ::EmProperty   ... three-photon transition strength S^(3) [atomic units]; NOT a cross
                                          section, see the header of this file.
"""
struct  Line_3pAbsorption
    initialLevel       ::Level
    finalLevel         ::Level
    omega1             ::Float64
    omega2             ::Float64
    omega3             ::Float64
    strength           ::EmProperty
end


# `Base.show(io::IO, line::MultiPhotonTransition.Line_3pAbsorption)`
#   ... prepares a proper printout of the variable line::MultiPhotonTransition.Line_3pAbsorption.
function Base.show(io::IO, line::MultiPhotonTransition.Line_3pAbsorption)
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "omega1:            $(line.omega1)  ")
    println(io, "omega2:            $(line.omega2)  ")
    println(io, "omega3:            $(line.omega3)  ")
    println(io, "strength:          $(line.strength)  ")
end


"""
`MultiPhotonTransition.allIntermediateLevels(mp::Multiplet)`  or
`MultiPhotonTransition.allIntermediateLevels(gChannels::Array{AtomicState.GreenChannel,1})`
    ... to return ALL intermediate levels of the given representation, irrespective of their symmetry; an
        Array{Level,1} is returned.

        The two-photon routines ask for the levels of ONE symmetry, because the coupled form fixes that symmetry
        in advance. The m-sum form below does not work with a predetermined symmetry: it lets the 3-j symbols and
        the parity rule decide which intermediate levels contribute, so it needs the whole set.
"""
function allIntermediateLevels(mp::Multiplet)
    return( mp.levels )
end

"""
`MultiPhotonTransition.allIntermediateLevels(gChannels::Array{AtomicState.GreenChannel,1})`
    ... collects the levels of ALL Green channels, of every symmetry. The unfiltered counterpart of
        `intermediateLevels`, needed by the three-photon amplitude, whose two intermediate sums do not share one
        total symmetry. A levels::Array{Level,1} is returned.
"""
function allIntermediateLevels(gChannels::Array{AtomicState.GreenChannel,1})
    levels = Level[]
    for  ch in gChannels    append!(levels, ch.gMultiplet.levels)    end
    return( levels )
end


"""
`MultiPhotonTransition.isStepAllowed(symA::LevelSymmetry, mp::EmMultipole, symB::LevelSymmetry)`
    ... to decide whether a one-photon step of multipolarity mp can connect the two level symmetries; a
        Bool is returned. Both the triangle rule and the parity rule are applied.

        THIS GUARD IS NOT OPTIONAL. `AngularMomentum.Wigner_3j` raises a DomainError for arguments that violate
        its own domain rather than returning zero, so forbidden steps must be filtered out before they are
        reached, not afterwards.
"""
function isStepAllowed(symA::LevelSymmetry, mp::EmMultipole, symB::LevelSymmetry)
    if  !AngularMomentum.isTriangle(symA.J, AngularJ64(mp.L), symB.J)                      return( false )   end
    if  !AngularMomentum.parityEmMultipolePi(symB.parity, mp, symA.parity)                 return( false )   end
    return( true )
end


"""
`MultiPhotonTransition.computeStrength_3pAbsorption(finalLevel::Level, initialLevel::Level,
                            omegas::Array{Float64,1}, gauge::EmGauge, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings; noOrderings::Int64=6)`
    ... to compute the three-photon transition strength S^(3) = 1/(2J_i+1) sum_M |A(M)|^2 for three parallel
        linearly-polarized beams of the given energies; a strength::Float64 is returned.

        THE AMPLITUDE, for one time ordering (a, b, c) of the three photons and with all polarizations along z,
        so that only the q = 0 components contribute and M is conserved:

            A(M) = sum_(nu1,nu2) (-1)^(J_f-M) 3j(J_f, L_c, J_nu2; -M, 0, M)
                                 (-1)^(J_nu2-M) 3j(J_nu2, L_b, J_nu1; -M, 0, M)
                                 (-1)^(J_nu1-M) 3j(J_nu1, L_a, J_i; -M, 0, M)
                                 <f||O(L_c,w_c)||nu2> <nu2||O(L_b,w_b)||nu1> <nu1||O(L_a,w_a)||i>
                                 / [(E_i + w_a - E_nu1) (E_i + w_a + w_b - E_nu2)]

        summed coherently over all 3! = 6 orderings. Note the two denominators: the intermediate state after the
        first absorption sits at E_i + w_a, the one after the second at E_i + w_a + w_b.

        `noOrderings` exists only for the acceptance test: setting it to 1 keeps the single ordering (1,2,3), so
        that the monochromatic limit can be checked to give exactly 6 times that amplitude.
"""
function computeStrength_3pAbsorption(finalLevel::Level, initialLevel::Level, omegas::Array{Float64,1},
                                      gauge::EmGauge, grid::Radial.Grid, settings::MultiPhotonTransition.Settings;
                                      noOrderings::Int64=6)
    nuLevels = MultiPhotonTransition.allIntermediateLevels(settings.intermediateStates)
    symi     = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symf     = LevelSymmetry(finalLevel.J,   finalLevel.parity)
    ji       = Basics.twice(initialLevel.J) / 2;    jf = Basics.twice(finalLevel.J) / 2
    # the reduced matrix elements are re-used across orderings and are cached on (levels, omega-index, multipole)
    cache    = Dict{Tuple{Int64,Int64,Int64,EmMultipole},ComplexF64}()
    function redME(levA::Level, levB::Level, iw::Int64, mp::EmMultipole)
        key = (levA.index, levB.index, iw, mp)
        if  haskey(cache, key)    return( cache[key] )    end
        wa  = PhotoEmission.amplitude(Absorption(), mp, gauge, omegas[iw], levA, levB, grid,
                                      display=false, printout=false)
        cache[key] = wa
        return( wa )
    end
    # M is conserved through the chain, so a single list serves; the 3-j symbols vanish outside it anyway
    Mlist    = collect(-min(ji,jf) : 1. : min(ji,jf))
    amps     = zeros(ComplexF64, length(Mlist))
    orders   = [(1,2,3), (1,3,2), (2,1,3), (2,3,1), (3,1,2), (3,2,1)][1:noOrderings]

    for  (ia, ib, ic) in orders
        wa = omegas[ia];    wb = omegas[ib];    wc = omegas[ic]
        for  mpa in settings.multipoles
            for  mpb in settings.multipoles
                for  mpc in settings.multipoles
                    for  nu1 in nuLevels
                        sym1 = LevelSymmetry(nu1.J, nu1.parity)
                        if  !MultiPhotonTransition.isStepAllowed(sym1, mpa, symi)    continue    end
                        denom1 = initialLevel.energy + wa - nu1.energy
                        if  abs(denom1) < settings.selfTolerance    continue    end
                        j1 = Basics.twice(nu1.J) / 2
                        for  nu2 in nuLevels
                            sym2 = LevelSymmetry(nu2.J, nu2.parity)
                            if  !MultiPhotonTransition.isStepAllowed(sym2, mpb, sym1)    continue    end
                            if  !MultiPhotonTransition.isStepAllowed(symf, mpc, sym2)    continue    end
                            denom2 = initialLevel.energy + wa + wb - nu2.energy
                            if  abs(denom2) < settings.selfTolerance    continue    end
                            j2 = Basics.twice(nu2.J) / 2
                            red = redME(finalLevel, nu2, ic, mpc) * redME(nu2, nu1, ib, mpb) *
                                  redME(nu1, initialLevel, ia, mpa) / (denom1 * denom2)
                            if  red == 0.    continue    end
                            for  (iM, M) in enumerate(Mlist)
                                if  abs(M) > j1  ||  abs(M) > j2    continue    end
                                wx = (-1.0)^(jf - M) * convert(Float64, AngularMomentum.Wigner_3j(jf, mpc.L, j2, -M, 0, M)) *
                                     (-1.0)^(j2 - M) * convert(Float64, AngularMomentum.Wigner_3j(j2, mpb.L, j1, -M, 0, M)) *
                                     (-1.0)^(j1 - M) * convert(Float64, AngularMomentum.Wigner_3j(j1, mpa.L, ji, -M, 0, M))
                                amps[iM] = amps[iM] + wx * red
                            end
                        end
                    end
                end
            end
        end
    end

    strength = 0.
    for  amp in amps    strength = strength + abs(amp)^2    end
    strength = strength / (Basics.twice(initialLevel.J) + 1)

    return( strength )
end


"""
`MultiPhotonTransition.determineOmegas_3pAbsorption(energy::Float64, scheme::ThreePhotonAbsorptionScheme)`
    ... to determine the three photon energies (in atomic units) from the scheme and the transition energy; an
        omegas::Array{Float64,1} is returned.

        omega1 = omega2 = 0 selects the MONOCHROMATIC case, all three photons carrying energy/3. Otherwise the
        two given energies are read in the user-selected units, like every other photon energy in JAC, and the
        third follows from energy conservation.
"""
function determineOmegas_3pAbsorption(energy::Float64, scheme::ThreePhotonAbsorptionScheme)
    if  scheme.omega1 == 0.  &&  scheme.omega2 == 0.
        return( [energy/3, energy/3, energy/3] )
    end
    omega1 = Defaults.convertUnits("energy: to atomic", scheme.omega1)
    omega2 = Defaults.convertUnits("energy: to atomic", scheme.omega2)
    omega3 = energy - omega1 - omega2
    if  omega1 <= 0.  ||  omega2 <= 0.  ||  omega3 <= 0.
        error("\n\nMultiPhotonTransition: the three photon energies must all be positive, but omega1 = $omega1, " *
              "omega2 = $omega2 and omega3 = $omega3 a.u. for a transition energy of $energy a.u.\n"              *
              ">>> Note that omega1 = omega2 = 0 selects the monochromatic case.\n")
    end
    return( [omega1, omega2, omega3] )
end


"""
`MultiPhotonTransition.computeLines(scheme::ThreePhotonAbsorptionScheme, finalMultiplet::Multiplet,
                            initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings; output=true)`
    ... to compute the three-photon absorption strengths for three parallel linearly-polarized beams; a list of
        lines::Array{MultiPhotonTransition.Line_3pAbsorption,1} is returned.
"""
function  computeLines(scheme::ThreePhotonAbsorptionScheme, finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                       grid::Radial.Grid, settings::MultiPhotonTransition.Settings; output=true)
    println("")
    printstyled("MultiPhotonTransition.computeLines(::ThreePhotonAbsorptionScheme): The computation of amplitudes starts now ... \n",
                color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    println("")
    newLines = MultiPhotonTransition.Line_3pAbsorption[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)    continue    end
            energy = fLevel.energy - iLevel.energy + settings.photonEnergyShift
            if  energy <= 0.    continue    end
            omegas = MultiPhotonTransition.determineOmegas_3pAbsorption(energy, scheme)
            if  Basics.UseCoulomb  in  settings.gauges
                    sCou = MultiPhotonTransition.computeStrength_3pAbsorption(fLevel, iLevel, omegas,
                                                        EmGauge("Coulomb"), grid, settings)
            else    sCou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    sBab = MultiPhotonTransition.computeStrength_3pAbsorption(fLevel, iLevel, omegas,
                                                        EmGauge("Babushkin"), grid, settings)
            else    sBab = 0.
            end
            push!( newLines, MultiPhotonTransition.Line_3pAbsorption(iLevel, fLevel, omegas[1], omegas[2],
                                                                     omegas[3], EmProperty(sCou, sBab)) )
        end
    end
    # Print all results to screen
    MultiPhotonTransition.displayResults_3pAbsorption(stdout, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    MultiPhotonTransition.displayResults_3pAbsorption(iostream, newLines)     end
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`MultiPhotonTransition.displayResults_3pAbsorption(stream::IO,
                            lines::Array{MultiPhotonTransition.Line_3pAbsorption,1})`
    ... to display the three-photon absorption strengths of the selected lines. A neat table is printed but
        nothing is returned otherwise.
"""
function  displayResults_3pAbsorption(stream::IO, lines::Array{MultiPhotonTransition.Line_3pAbsorption,1})
    nx = 118
    println(stream, " ")
    println(stream, "  Three-photon absorption by three parallel linearly-polarized beams:")
    println(stream, " ")
    println(stream, "  The three-photon transition STRENGTH  S^(3) = 1/(2J_i+1) sum_M |A(M)|^2  is given in atomic units.")
    println(stream, "  It is NOT a generalized cross section: that would need a normalisation of order F^3, and the")
    println(stream, "  two-photon absorption normalisation of this module has never been derived either.")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(10, "omega1"; na=3);   sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "omega2"; na=3);   sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "omega3"; na=4);   sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(30, "Cou --   S^(3) [a.u.]   -- Bab"; na=3);     sb = sb * TableStrings.hBlank(33)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa   = "";   isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                     fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega1))  * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega2))  * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega3))  * "      "
        sa = sa * @sprintf("%.5e", line.strength.Coulomb)    * "      "
        sa = sa * @sprintf("%.5e", line.strength.Babushkin)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">>> A strength of exactly 0 is a SELECTION RULE, not a failure: three E1 photons change the " *
                    "parity, so a\n>>> transition between levels of the SAME parity cannot be driven by them.")
    return( nothing )
end
