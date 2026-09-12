#
# probe-eolOrbitalEnergy.jl   --   groundwork for priority item 11.
#
# Every orbital of an EOL or RAS basis carries `energy = 0.0` exactly, because a rotation-optimized orbital is
# not the eigenfunction of any one-particle operator.  Item 11 asks for a DEFINED energy, names the diagonal
# Lagrange multiplier eps_a = <a|F|a> as the candidate, and warns that the one-particle expectation <a|h_D|a>
# is NOT the orbital energy and is the plausible-looking wrong answer.
#
# THIS PROBE DECIDES BETWEEN THEM WITHOUT TOUCHING THE SOLVER.  It builds both quantities from the orbitals
# alone and compares them with the AL orbital energies of the same closed-shell system -- the check item 11
# itself proposes: "the AL and EOL orbital energies of a closed-shell case should agree to the size of the
# correlation they differ by, not by hundreds of Hartree".
#
# WHAT IT MEASURES, AND WHAT IT DOES NOT.  Three columns are printed: the AL orbital energy (the reference any
# definition has to reproduce to within a correlation-sized amount), the 0.0 the EOL solver stores today, and
#
#   <a|h_D|a>            the one-particle (Dirac + nucleus) expectation alone.
#
# THE FOCK EXPECTATION <a|F|a> IS NOT COMPUTED HERE.  Building it needs the tensor caches that the rotation
# solver does not construct, which is the implementation work item 11 is actually about; this probe exists to
# settle the CHEAPER question first -- whether the one-particle expectation could serve instead -- because that
# route needs nothing new and would be the tempting shortcut.  It cannot serve: see the numbers it prints.
# So this probe EXCLUDES an alternative; it does not yet CONFIRM <a|F|a>.
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
    @printf("%-10s %16s %16s %16s\n", "subshell", "AL energy", "EOL stored", "<a|h_D|a>")
    for sh in bAL.subshells
        eAL = bAL.orbitals[sh].energy
        eST = haskey(bEOL.orbitals, sh) ? bEOL.orbitals[sh].energy : NaN
        h1  = haskey(bEOL.orbitals, sh) ? oneParticle(bEOL.orbitals[sh], grid, nucPot) : NaN
        @printf("  %-8s %16.5f %16.5f %16.5f\n", string(sh), eAL, eST, h1)
    end
    println("\n  READ THE THIRD COLUMN AGAINST THE FIRST: <a|h_D|a> omits ALL electron-electron interaction, so it")
    println("  must come out FAR BELOW the true orbital energy -- that is exactly why item 11 warns against it.")
    println("  The second column is the 0.0 the solver stores today.")
end
