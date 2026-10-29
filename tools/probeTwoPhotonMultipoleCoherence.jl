#
#  [S22]  Does the two-photon emission prefactor and its CHANNEL GROUPING generalise to K != 0 and to
#         mixed multipoles?  Challenge S22.  THIS PROBE CHANGES NOTHING IN `src/`: it re-forms the very
#         same `wk` that `MultiPhotonTransition.computeEnergyDiffCs` forms, out of the reduced channels
#         that `Line_2pEmission` already carries, and then asks two questions the production code cannot
#         answer about itself.
#
#  WHAT IT FOUND, 29-Oct-2026 -- and the second half outgrew the question it was written for.
#
#  RUN A -- [S22-ODDK].  ONE QUESTION: does an ODD K vanish at EVERY energy sharing, or only at the
#  symmetric point omega1 = omega2?  The challenge's own record states the strong form ("odd K vanishes by
#  Bose exchange symmetry whenever BOTH photons carry the SAME multipole type, for ANY transition"), and
#  the suspicion was that this had been over-generalised from the midpoint.  MEASURED, H-like 3d_3/2 -> 1s
#  in pure E1E1 over ten sharings:  K = 1 sits at 1e-21 ... 1e-15 against a K = 2 of order 1e-12 at EVERY
#  sharing, not merely at the midpoint.  SO THE SUSPICION WAS WRONG AND THE RECORD STANDS;  the practical
#  conclusion it draws -- that a pure-E1E1 case cannot exercise odd K at all -- is confirmed.
#
#  RUN B -- [S22-CHAN].  ONE QUESTION: of the three M2-involving groups that came back EXACTLY zero in the
#  documented dead end, which are zero because the PHYSICS is zero and which because the code never finds
#  them?  ANSWER: ALL THREE WERE CODE, in three different ways, and fixing them is what this run now
#  demonstrates rather than diagnoses.  Before (H 2s -> 1s, multipoles [E1, M2], sharing 1):
#
#      group        stored     max|amplitude|   retrieved by the production lookup
#      (E1,E1)         4        1.20e-06        found
#      (E1,M2)         2        0.00e+00        found, and exactly zero
#      (M2,E1)         2        0.00e+00        found, and exactly zero
#      (M2,M2)         8        1.23e-17        NOT FOUND -- warning, returns 0
#
#  The three causes, none of which is a selection rule:
#    (1) `InteractionStrength.MabEmission` branched on the GAUGE ALONE and never looked at the multipole, so
#        the M2 factor of a mixed pair, filed under an electric gauge, was computed with the ELECTRIC E2
#        formula.  E2 and M2 carry OPPOSITE parity, so the substituted operator was forbidden for exactly
#        the orbital pairs the right one allows -- hence exactly zero rather than merely wrong.  Proof:
#        M2/Coulomb equalled E2/Coulomb to all 17 digits, and M1/Magnetic = +3.8150591874e-12 where
#        M1/Coulomb = 0, which calibrates the detector.
#    (2) `determineSharingsAndChannels` filed a pure-magnetic pair under the `Magnetic` gauge key, which
#        `getReducedAmplitudeEmission` never asks for -- its one-line `Magnetic` disjunct had been commented
#        out rather than finished.  Built, filled, then reported missing.
#    (3) The same pair was filed once per requested gauge, i.e. stored twice over (8 channels against 4).
#
#  AFTER (same case, same sharing) -- every group live, correctly keyed, and no duplicates:
#
#      gauge key    mp1  mp2   channels   max|amplitude|
#      Babushkin    E1   E1         4     3.037905e-06
#      Babushkin    E1   M2         2     1.648146e-11
#      Babushkin    M2   E1         2     2.093424e-12
#      Babushkin    M2   M2         4     1.759480e-17
#      Coulomb      E1   E1         4     9.516940e-07
#      Coulomb      E1   M2         2     3.798879e-12
#      Coulomb      M2   E1         2     2.298035e-12
#      Coulomb      M2   M2         4     1.759480e-17
#
#  READ THOSE MAGNITUDES AS A STRUCTURAL RESULT AND NOT AS PHYSICS.  This case's intermediate basis is np
#  only, so it does not span the j = 5/2 symmetries an M2 factor can reach -- the run says so itself, with a
#  "the intermediate basis does not span this symmetry" warning per missing path.  What is demonstrated is
#  that the channels EXIST, are keyed where the lookup reaches them, and are no longer duplicated;  a
#  mixed-multipole RATE needs nd/nf intermediates and a system where E1M2 is not a 1e-5 relativistic
#  afterthought.  The literature anchor for that is Derevianko & Johnson, PRA 56, 1288 (1997): He-like
#  2^3S_1 -> 1^1S_0, where E1M2 is 18 % of the 2E1 rate at Z = 30.
#
#  SO S22'S OWN QUESTION -- whether the prefactor's grouping is coherent or incoherent ACROSS (mp1, mp2) --
#  IS NOW TESTABLE AND WAS NOT BEFORE, because until this was fixed JAC computed no mixed-multipole
#  two-photon amplitude at all.  That question is still open.
#
#  RUN:  julia --project=. tools/probeTwoPhotonMultipoleCoherence.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")

## `wk` EXACTLY AS THE PRODUCTION CODE FORMS IT: the Jsym sum is coherent and inside the modulus, and the
## two photon orderings are added with the exchange phase (-1)^(L1+L2-K).  Reproduced here rather than
## called, because `computeEnergyDiffCs` returns only the summed cross section.
function wkOf(line, sharing, K, mp1, mp2, gauge)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)
    wk   = Complex(0.)
    for  Jsym  in  AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
        wk = wk + MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, sharing.omega2, mp2, Jsym,
                                        sharing.omega1, mp1, line.initialLevel, gauge, line.sharings) +
                  (-1.0)^( mp1.L + mp2.L - Basics.twice(K)/2 ) *
                  MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, sharing.omega1, mp1, Jsym,
                                        sharing.omega2, mp2, line.initialLevel, gauge, line.sharings)
    end
    return( wk )
end

function runCase(tag, nm, grid, initConfs, finalConfs, interConfs, mps, noSharings, scf)
    asf      = AsfSettings(AsfSettings(); scField=scf, gridStopper=false)
    interRep = Representation("intermediate levels", nm, grid, interConfs, MeanFieldMultiplet(MeanFieldSettings(scf, false)))
    interMp  = redirect_stdout(devnull) do
        generate(interRep, output=true)[ResultKeys.MeanFieldMultiplet]
    end
    mpSet = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                             MultiPhotonTransition.AbstractMultiPhotonProperty[
                                 MultiPhotonTransition.EnergyDiffCs()], noSharings ),
                multipoles = mps, gauges = [UseCoulomb, UseBabushkin],
                intermediateStates = interMp, calcOverview = false, lineSelection = LineSelection() )
    wa = redirect_stdout(devnull) do
        perform( Atomic.Computation(Atomic.Computation(), name=tag, grid=grid, nuclearModel=nm,
                                   initialConfigs=initConfs, finalConfigs=finalConfs,
                                   initialAsfSettings=asf, finalAsfSettings=asf,
                                   processSettings=mpSet); output=true )
    end
    return( wa[MultiPhotonTransition.Settings] )
end

############################################################  RUN A -- [S22-ODDK]
println("\n", "="^118, "\n  RUN A  [S22-ODDK]   H-like 3d_3/2 -> 1s_1/2, PURE E1E1, resolved in the energy sharing.")
println("  ONE QUESTION: is odd K zero everywhere, or only at omega1 = omega2?\n", "="^118)

nmA   = Nuclear.Model(1.0, PointNucleus())
gridA = Basics.recommendedGrid(Dict(Shell("1s") => 1, Shell("3d") => 1, Shell("8p") => 1), 1.0)
setDefaults("standard grid", gridA)
interA = [Configuration("$(n)p") for n = 2:8]
linesA = runCase("S22-oddK", nmA, gridA, [Configuration("3d")], [Configuration("1s")], interA, [E1], 10, Basics.NuclearField())

for  line in linesA
    Basics.twice(line.initialLevel.J) == 3  ||  continue
    @printf("\n  line %s -> %s,  total omega = %.4f eV,  %d sharings\n",
            string(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)),
            string(LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)),
            (line.initialLevel.energy-line.finalLevel.energy)*27.211386245988, length(line.sharings))
    println("     omega1/omega_tot      |wk(K=1)|^2        |wk(K=2)|^2      odd-K share")
    for  sh  in  line.sharings
        tot = sh.omega1 + sh.omega2
        w1  = abs(wkOf(line, sh, AngularJ64(1), E1, E1, Basics.Coulomb))^2
        w2  = abs(wkOf(line, sh, AngularJ64(2), E1, E1, Basics.Coulomb))^2
        @printf("        %8.5f        %.6e      %.6e     %.3e\n", sh.omega1/tot, w1, w2, w1/(w1+w2))
    end
end

############################################################  RUN B -- [S22-CHAN]
println("\n\n", "="^118, "\n  RUN B  [S22-CHAN]   H 2s -> 1s with multipoles [E1, M2] -- the DOCUMENTED DEAD END, now instrumented.")
println("  ONE QUESTION: is every (mp1, mp2) group STORED, keyed where the production lookup reaches it, and stored")
println("  exactly once?  Compare against the two tables in this file's header.\n", "="^118)

nmB   = Nuclear.Model(1.0, PointNucleus())
gridB = Basics.recommendedGrid(Dict(Shell("1s") => 1, Shell("2s") => 1, Shell("8p") => 1), 1.0)
setDefaults("standard grid", gridB)
interB = vcat([Configuration("$(n)p") for n = 2:8], [Configuration("$(n)d") for n = 3:6],
              [Configuration("$(n)f") for n = 4:6])
linesB = runCase("S22-chan", nmB, gridB, [Configuration("2s")], [Configuration("1s")], interB, [E1, M2], 4,
                 Basics.NuclearField())

for  line in linesB
    sh = line.sharings[1]
    @printf("\n  line, total omega = %.4f eV;  sharing 1 at omega1/omega_tot = %.5f;  %d stored channels\n",
            (line.initialLevel.energy-line.finalLevel.energy)*27.211386245988,
            sh.omega1/(sh.omega1+sh.omega2), length(sh.channels))

    ## WHAT IS ACTUALLY STORED, grouped by the gauge KEY it was filed under
    println("\n    stored channels, by (gauge key, mp1, mp2):   count   max|amplitude|")
    seen = Dict{Tuple{String,String,String},Vector{Float64}}()
    for ch in sh.channels
        k = (string(ch.gauge), string(ch.multipole1), string(ch.multipole2))
        push!( get!(seen, k, Float64[]), abs(ch.amplitude) )
    end
    for k in sort(collect(keys(seen)))
        @printf("      %-12s %-4s %-4s                     %3d     %.6e\n", k[1], k[2], k[3],
                length(seen[k]), maximum(seen[k]))
    end

    ## WHAT computeEnergyDiffCs ACTUALLY RETRIEVES for the same groups
    println("\n    what the production lookup returns, |wk|^2:")
    println("      mp1   mp2      K        Coulomb          Babushkin         Magnetic")
    for mp1 in [E1, M2],  mp2 in [E1, M2]
        for K in [AngularJ64(0), AngularJ64(1)]
            v = map(g -> abs(wkOf(line, sh, K, mp1, mp2, g))^2,
                    [Basics.Coulomb, Basics.Babushkin, Basics.Magnetic])
            @printf("      %-4s  %-4s   %-4s    %.6e    %.6e    %.6e\n",
                    string(mp1), string(mp2), string(K), v[1], v[2], v[3])
        end
    end
end
println()
