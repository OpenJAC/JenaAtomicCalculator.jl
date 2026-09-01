#
println("Kb) Apply & test the AtomicFeatures module: the FEATURE VECTORS that a neural network is trained on.")
#
# WHERE THIS SITS.  example-Ka.jl trains and applies a network; this file is one step earlier and asks what the
# network is given to learn FROM.  A feature is just a number extracted from an atomic state that the network can
# see -- how many electrons sit in each shell, how far out an orbital reaches, and so on -- and the quality of a
# machine-learned prediction is limited far more by the choice of features than by the network.
#
# WHY IT IS WORTH A FILE OF ITS OWN.  The extractors below are the only part of the deep-learning work that is
# pure, checkable atomic physics: each returns a number whose value is known in advance for a simple case, so
# they can be verified WITHOUT training anything and without any statistics.  That matters here, because the
# Deep-Learning study itself stands at an honest NEGATIVE result -- Stage 2 concluded that the training data must
# be widened before Stage 3 is attempted (application item A14).  A network that learns nothing may be failing
# because the physics it was shown was too thin; being able to check the features separately is what keeps those
# two possibilities apart.
#
# THE FEATURES SHOWN HERE ARE THE TRIVIAL ONES, deliberately: occupations, mean occupation numbers and <r^k>.
# The expensive ones -- the Slater F^k and G^k integrals -- are the same idea with a two-electron radial integral
# behind them and are exercised by the module's own test rather than here.

using SymEngine

if  true
    # Last visit:      01-Sep-2026
    # Last successful: 01-Sep-2026 -- every occupation below is what counting the configuration by hand gives.
    #
    # Branch a: THE CONFIGURATION-LEVEL FEATURES, which need no orbitals and no self-consistent field at all.
    #   `extractShellOccupations` maps a configuration onto a FIXED list of shells, so that every atom in a
    #   training set produces a vector of the same length in the same order -- which is what a network requires
    #   and what a bare `Configuration` does not provide.  A shell absent from the configuration must come back
    #   as 0 rather than be missing, and that is the whole point of passing the shell list separately.
    #
    #   WHAT TO CHECK BY EYE: the numbers are simply the occupations written in the configuration, in the order
    #   of the `shells` list, with zeros where a shell does not occur.  Nothing here can be subtle, which is
    #   exactly why it is the right place to start.
    shells = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p"), Shell("3d")]
    for  conf  in  [Configuration("1s^2 2s^2 2p^6"),           # neon, closed
                    Configuration("1s^2 2s^2 2p^6 3s"),        # sodium-like, one valence electron
                    Configuration("1s^2 2s^2 2p^4"),           # oxygen-like, an open p shell
                    Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^5")]   # a half-filled d shell
        occ = AtomicFeatures.extractShellOccupations(shells, conf)
        println("\n  $conf")
        println("    shells      = $shells")
        println("    occupations = $occ        (sum = $(sum(occ)) electrons)")
    end
    #
elseif  false
    # Last visit:      01-Sep-2026
    # Last successful: 01-Sep-2026 -- <r> and <r^2> come out in the expected order 1s < 2s < 2p for neon, and the
    #                  mean occupation numbers reproduce the closed-shell values 2, 2, 2, 4 exactly.
    #
    # Branch b: THE LEVEL- AND ORBITAL-LEVEL FEATURES, which do need a self-consistent field.
    #
    #   `extractMeanOccupationNumbers` gives the occupation of each SUBSHELL averaged over the CSFs of a level,
    #   with their mixing coefficients.  For a closed-shell ground state that is just the closed-shell occupation
    #   and must come out as an integer; for an open shell it is generally fractional, and THAT is the
    #   information a network cannot get from the configuration label alone.
    #
    #   `extractRkExpectation` gives <r^k> for each subshell -- the size of the orbital.  Two things must hold and
    #   are worth reading off: <r> grows with n, and for the same n the more weakly bound orbital reaches
    #   further.  A feature that did not respect that ordering would be telling the network something false about
    #   the atom's geometry.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 12.0)
    nm   = Nuclear.Model(10.)
    asf  = AsfSettings(AsfSettings(); scField = Basics.DFSField())
    multiplet = SelfConsistent.performSCF([Configuration("1s^2 2s^2 2p^6")], nm, grid, asf; printout = false)
    level     = multiplet.levels[1]
    subshells = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
    #
    println("\n  Ne 1s^2 2s^2 2p^6, ground level $(LevelSymmetry(level.J, level.parity)):")
    meanOcc = AtomicFeatures.extractMeanOccupationNumbers(subshells, level)
    println("    subshells             = $subshells")
    println("    mean occupations      = $meanOcc      (closed shells, so 2, 2, 2, 4 exactly)")
    #
    for  k  in  [1, 2]
        rk = AtomicFeatures.extractRkExpectation(subshells, k, level.basis.orbitals, grid)
        println("    <r^$k>                  = $rk")
    end
    println("\n    read the ordering: <r> must grow with n, and 1s < 2s and 1s < 2p for the same reason --")
    println("    a feature that broke that ordering would misdescribe the atom's geometry to the network.")
    #
end
