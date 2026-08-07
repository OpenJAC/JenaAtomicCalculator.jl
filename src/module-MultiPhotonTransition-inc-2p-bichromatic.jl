#
# Two-photon absorption by BI-CHROMATIC photons, i.e. two beams of different frequency.
#
# NEW FILE, 06-Aug-2026. `TwoPhotonAbsorptionBichromatic` had existed as a process type since the module was
# written, was offered in the docstring of the abstract type, and had its own default constructor -- but no
# `-inc-` file and therefore no `computeLines` method at all, so selecting it died in a MethodError with
# nothing to say why. This file closes that gap.
#
# THE PHYSICS DIFFERENCE from the monochromatic case is one line, and it is the reason this case is worth
# having beyond its own applications: the two photons come from distinguishable beams, so
#
#     omega1 = scheme.omegaLess          omega2 = (E_f - E_i) - omega1
#
# are FIXED and unequal, the two orderings of the photons are physically distinct, and the rate is
# W = sigma^(2) * F_1 * F_2 with no combinatorial factor to argue about. The single-beam case, by contrast,
# involves indistinguishable photons and therefore a convention. Requiring the bichromatic result to go over
# into the monochromatic one as omega1 -> omega2 is what FIXES that convention rather than leaving it to be
# guessed -- see the units note in MultiPhotonTransition.Settings.
#

"""
`MultiPhotonTransition.computeLines(scheme::TwoPhotonAbsorptionBichromaticScheme, finalMultiplet::Multiplet,
                            initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings; output=true)`
    ... to compute the two-photon absorption lines for bi-chromatic photons; a list of
        lines::Array{MultiPhotonTransition.Line_2pAbsorptionMonochromatic,1} is returned.

        STAGE 1 (06-Aug-2026): the line/channel machinery and the reduced amplitudes are shared with the
        monochromatic case, which already carries two independent multipoles mp1, mp2 and their two orderings;
        what differs is only that the two photon energies are unequal. Since the stored
        `Line_2pAbsorptionMonochromatic` carries a single `omega`, the amplitudes below are evaluated at
        `omegaLess` and the asymmetric-energy pieces are NOT yet formed -- this method therefore raises an
        error rather than returning a half-correct number. The scaffolding is in place; the physics is Stage 2.
"""
function  computeLines(scheme::TwoPhotonAbsorptionBichromaticScheme, finalMultiplet::Multiplet,
                       initialMultiplet::Multiplet, grid::Radial.Grid,
                       settings::MultiPhotonTransition.Settings; output=true)
    error("\n\nMultiPhotonTransition.computeLines():  STOP -- TwoPhotonAbsorptionBichromaticScheme is scaffolded "  *
          "but not yet complete.\n"                                                                                *
          ">>> What is still missing:\n"                                                                           *
          "    (1) a Line/Channel type carrying TWO photon energies. Line_2pAbsorptionMonochromatic stores one\n"   *
          "        `omega`, which is sufficient only while both photons share it;\n"                               *
          "    (2) the reduced amplitude evaluated with the two DIFFERENT denominators\n"                          *
          "        (E_i + omega1 - E_nu) and (E_i + omega2 - E_nu), one per ordering of the two photons;\n"         *
          "    (3) the cross section as W = sigma^(2) * F_1 * F_2 rather than F^2.\n"                               *
          ">>> Note that (3) is the UNAMBIGUOUS normalisation and is meant to fix the single-beam convention of\n" *
          "    TwoPhotonAbsorptionScheme by the limit omega1 -> omega2; see the units note in Settings.\n"          *
          ">>> Use TwoPhotonAbsorptionScheme() meanwhile.\n")
end
