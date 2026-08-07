#
# Three-photon excitation and decay between bound levels.
#
# NOT YET IMPLEMENTED. This file exists so that the scheme hierarchy is complete and so that selecting a
# three-photon scheme fails with a message that says what is missing, rather than dying in a MethodError far
# from the cause -- which is what happened with TwoPhotonAbsorptionBichromatic, whose scheme type existed while
# its computeLines did not.
#

"""
`MultiPhotonTransition.computeLines(scheme::ThreePhotonEmissionScheme, finalMultiplet::Multiplet,
                            initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings; output=true)`
    ... to compute three-photon emission lines. NOT YET IMPLEMENTED; an informative error is raised.
"""
function  computeLines(scheme::ThreePhotonEmissionScheme, finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                       grid::Radial.Grid, settings::MultiPhotonTransition.Settings; output=true)
    error("\n\nMultiPhotonTransition.computeLines():  STOP -- ThreePhotonEmissionScheme is not yet implemented.\n" *
          ">>> What is missing, in the order it would have to be built:\n"                                        *
          "    (1) a THIRD-order amplitude, i.e. a double sum over two sets of intermediate states\n"             *
          "        <f|O(mp3)|nu2> <nu2|O(mp2)|nu1> <nu1|O(mp1)|i> / [(E_i + w1 - E_nu1)(E_i + w1 + w2 - E_nu2)]\n"*
          "        with all 3! = 6 orderings of the photons, against 2! = 2 for the two-photon case;\n"           *
          "    (2) energy sharings on a two-dimensional simplex (w1 + w2 + w3 = E_i - E_f) rather than the\n"     *
          "        one-dimensional Gauss-Legendre line used for two photons;\n"                                   *
          "    (3) the coupling of THREE photon multipoles, so the rank K is no longer fixed by\n"                *
          "        oplus(J_f, J_i) alone and an intermediate coupling rank must be summed over.\n"                *
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
          ">>> Use TwoPhotonAbsorptionScheme() or TwoPhotonAbsorptionBichromaticScheme() meanwhile.\n")
end
