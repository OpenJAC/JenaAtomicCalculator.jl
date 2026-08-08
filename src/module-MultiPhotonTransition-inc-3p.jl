#
# Three-photon excitation and decay between bound levels.
#
# THE AMPLITUDE IS NOT IMPLEMENTED, and this file does not pretend otherwise: computeLines still refuses, with a
# message naming what would have to be built. What IS implemented, 08-Aug-2026, is the piece that can be built
# and VERIFIED on its own -- the ENERGY SHARINGS ON THE TWO-DIMENSIONAL SIMPLEX.
#
# WHY THAT PIECE FIRST. For two photons the sharings are a Gauss-Legendre line: omega1 in [0, E], omega2 fixed by
# energy conservation. For three photons the domain is the triangle
#
#     omega1 + omega2 + omega3 = E,      omega1, omega2, omega3 > 0
#
# and the rate is differential in TWO of them. Getting that quadrature right is a self-contained problem with an
# EXACT answer to check against, which is what makes it worth doing before the third-order amplitude: the moments
# of a simplex are known in closed form,
#
#     Int omega1^a omega2^b omega3^c domega1 domega2  =  E^(a+b+c+2) * a! b! c! / (a+b+c+2)!
#
# so a = b = c = 0 must give E^2/2, (1,0,0) must give E^3/6, and (1,1,1) must give E^5/120. Three sum rules of
# increasing polynomial degree, all parameter-free. If the weights are wrong, the eventual three-photon total
# rate is wrong by a factor that no amount of checking the amplitude would ever reveal -- exactly the class of
# silent normalisation error that has cost this module time before.
#
# THE TRANSFORM is the standard collapsed-coordinate (Duffy) map of the square onto the triangle,
#
#     omega1 = E u,     omega2 = E (1-u) v,     omega3 = E (1-u)(1-v),     u, v in [0,1]
#
# with Jacobian domega1 domega2 = E^2 (1-u) du dv. A tensor Gauss-Legendre rule in (u,v) then integrates any
# polynomial exactly up to the order of the rule, including the (1-u) weight. `noSharings` is the number of
# points PER DIRECTION, so the grid carries noSharings^2 sharings -- stated because it is the one place a user
# could reasonably expect it to mean the total.
#
# ONE HONEST LIMITATION, stated rather than discovered later: the Duffy grid is EXACT but not PERMUTATION
# SYMMETRIC -- omega1 plays a distinguished role in the transform. For the integrated rate that is irrelevant
# (exactness is exactness), but a three-photon spectrum plotted on these nodes will not LOOK symmetric under
# exchanging the photons, even though the underlying function is. The symmetry check that served the two-photon
# case so well (blocker A1 was found through it) therefore has to be applied to the function, evaluated at
# permuted points, and not to the tabulated node values.
#


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
            ## the Jacobian domega1 domega2 = E^2 (1-u) du dv
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
    #
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
    #
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
    ## calcOverview shows what exists rather than refusing outright. The sharings need nothing but the transition
    ## energy, so they can be inspected and checked long before any amplitude is written -- which is the point of
    ## having built them first.
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


"""
`MultiPhotonTransition.computeLines(scheme::ThreePhotonAbsorptionScheme, finalMultiplet::Multiplet,
                            initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings; output=true)`
    ... to compute three-photon absorption lines. NOT YET IMPLEMENTED; an informative error is raised.
"""
function  computeLines(scheme::ThreePhotonAbsorptionScheme, finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                       grid::Radial.Grid, settings::MultiPhotonTransition.Settings; output=true)
    error("\n\nMultiPhotonTransition.computeLines():  STOP -- ThreePhotonAbsorptionScheme is not yet implemented.\n" *
          ">>> It needs the same third-order amplitude as ThreePhotonEmissionScheme (see that error message),\n"     *
          "    plus a generalized cross section of order F^3, i.e. units cm^6 s^2 rather than cm^4 s.\n"            *
          ">>> Note that three photons of FIXED energies do not need sharings at all -- the sharings belong to\n"   *
          "    the emission case, where only the TOTAL energy is fixed.\n"                                          *
          ">>> Use TwoPhotonAbsorptionScheme() or TwoPhotonAbsorptionBichromaticScheme() meanwhile.\n")
end
