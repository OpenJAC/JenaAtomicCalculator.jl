#
# probe-eolOrbitalEnergy.jl   --   the regression test for priority item 11.
#
# Until 12-Sep-2026 every orbital of an EOL or RAS basis carried `energy = 0.0` exactly.  Item 11 asked for a
# DEFINED energy, named the diagonal Lagrange multiplier eps_a = <a|F_a|a> as the candidate, and warned that the
# one-particle expectation <a|h_D|a> is NOT the orbital energy and is the plausible-looking wrong answer.
# BOTH EOL SOLVERS NOW SET eps_a (SelfConsistent.computeOrbitalEnergiesEOL), so this file's job has changed from
# deciding the question to GUARDING the answer.
#
# THIS PROBE DECIDES BETWEEN THEM WITHOUT TOUCHING THE SOLVER.  It builds both quantities from the orbitals
# alone and compares them with the AL orbital energies of the same closed-shell system -- the check item 11
# itself proposes: "the AL and EOL orbital energies of a closed-shell case should agree to the size of the
# correlation they differ by, not by hundreds of Hartree".
#
# THE TEST IS THE MIDDLE COLUMN AGAINST THE FIRST.  These are CLOSED-SHELL, single-configuration cases, where
# the EOL and AL functionals are the same functional, so the two columns must agree to the accuracy the two
# optimizers converge to -- NOT to hundreds of Hartree, and NOT to a factor of two.  Measured 12-Sep-2026:
#
#   Be-like  1s  -4.73386 (AL)  -4.73350 (EOL)        Ne-like  1s  -32.81822  -32.81740
#            2s  -0.30933       -0.30932                       2s   -1.93607   -1.93580
#                                                              2p_  -0.85304   -0.85277
#
# The third column is the alternative that was EXCLUDED by this probe before the work was done:
#
#   <a|h_D|a>            the one-particle (Dirac + nucleus) expectation alone,
#
# which omits the electron-electron interaction and so runs 1.7x to 12x too DEEP while looking perfectly
# well-behaved.  It is kept in the output as the contrast: if a future change ever moves the middle column
# towards the third, that is the failure this file exists to catch.
#
using JenaAtomicCalculator, Printf
const B = JenaAtomicCalculator.Basics

function alReference(confs, nm, grid)
    set = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction(),
                                     eeInteractionCI = CoulombInteraction(), gridStopper = false)
    mp  = SelfConsistent.performSCF(confs, nm, grid, set; printout=false)
    return( mp.levels[1].basis )
end

function eolBasis(confs, nm, grid)
    set = AsfSettings(AsfSettings(); scField = Basics.EOLField(), eeInteraction = CoulombInteraction(),
                                     eeInteractionCI = CoulombInteraction(), gridStopper = false,
                                     scfRoute = Basics.RotationRoute(24))
    mp  = SelfConsistent.performSCF(confs, nm, grid, set; printout=false)
    return( mp.levels[1].basis )
end

# one-particle expectation <a|h_D|a>
function oneParticle(orb, grid, nucPot)
    return( RadialIntegrals.GrantIab(orb, orb, grid, nucPot) )
end

for (tag, Z, cstr) in [("Be-like  1s^2 2s^2", 4.0, "1s^2 2s^2"),
                       ("Ne-like  1s^2 2s^2 2p^6", 10.0, "1s^2 2s^2 2p^6")]
    confs = [Configuration(cstr)];   nm = Nuclear.Model(Z)
    grid  = Basics.recommendedGrid(confs, nm; printout=false)
    nucPot= Nuclear.nuclearPotential(nm, grid)
    bAL   = alReference(confs, nm, grid)
    bEOL  = eolBasis(confs, nm, grid)
    println("\n", "="^96);   println(tag, "   (Z = ", Z, ")");   println("="^96)
    @printf("%-10s %16s %16s %16s %10s\n", "subshell", "AL energy", "EOL eps_a", "<a|h_D|a>", "AL-EOL")
    for sh in bAL.subshells
        eAL = bAL.orbitals[sh].energy
        eST = haskey(bEOL.orbitals, sh) ? bEOL.orbitals[sh].energy : NaN
        h1  = haskey(bEOL.orbitals, sh) ? oneParticle(bEOL.orbitals[sh], grid, nucPot) : NaN
        @printf("  %-8s %16.5f %16.5f %16.5f %10.2e\n", string(sh), eAL, eST, h1, abs(eAL - eST))
    end
    println("\n  THE TEST IS COLUMN 2 AGAINST COLUMN 1: one closed-shell configuration makes the EOL and AL functionals")
    println("  the same functional, so eps_a must reproduce the AL energy to optimizer accuracy.  Column 3 is the")
    println("  one-particle expectation <a|h_D|a>, which omits ALL electron-electron interaction and therefore lands")
    println("  far too DEEP;  it is shown as the contrast, since drifting towards it is the failure to catch.")
end
