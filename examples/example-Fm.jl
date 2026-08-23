
println("Fm) Cascade.DielectronicCaptureScheme: the resonant capture of a free electron into doubly-excited")
println("    levels, for the KLL group of helium-like carbon.")
println("    This is the FIRST of the twelve schemes of Table 1 of the 2024 EPJD cascade paper to get an")
println("    example file of its own.  It had none until now, and was reachable only indirectly, through the")
println("    resonant channel of Cascade.ElectronExcitationScheme.")

using Printf

grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 8.0)
setDefaults("standard grid", grid)
## Radial.Grid(true) must NOT be used: it carries hp = 0 and the continuum orbitals of an Auger step then
## cannot be generated at all.

captureScheme() = Cascade.DielectronicCaptureScheme(500.0, 0., 1, [Shell("1s")], [Shell("2s"), Shell("2p")],
                                                    [Shell("2s"), Shell("2p")])

if  true
    # Last successful:  23-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest
    #    had been resolved under 1.12.6 on the maintainer machine.  Re-run there to confirm.]
    #   1.3 s WARM; 9 Auger steps, 16 capture lines.
    #
    # Branch a: THE COMPUTATION.  A free electron is captured by helium-like carbon while a 1s electron is
    #   excited to n = 2, giving the 1s 2l 2l' doubly-excited levels of lithium-like carbon -- the KLL group.
    #   The scheme computes ONLY the capture, i.e. the Auger widths of those resonances by detailed balance;
    #   the radiative stabilization that would turn this into dielectronic recombination belongs to
    #   Cascade.DielectronicRecombinationScheme and is not computed here.
    #
    #   ON THE COST, since estimating it in advance is the discipline: 1.3 s WARM.  The first cascade of a
    #   Julia session pays the compilation of the whole path and took 30 s for exactly this run -- a factor of
    #   23.  Any timing quoted in this series is the warm one, and a cold first measurement is not a
    #   measurement of the physics.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fm-capture.sum")
    #
    wa = Cascade.Computation(Cascade.Computation(); name="KLL dielectronic capture of He-like C",
                             nuclearModel=Nuclear.Model(6.), grid=grid, approach=Cascade.AverageSCA(),
                             scheme=captureScheme(), initialConfigs=[Configuration("1s^2")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    println("\n  >> capture lines: ", length(wb["dielectronic-capture lines:"]))
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  23-Aug-2026
    #   [PROVENANCE: as branch a.]  ~2.2 s warm for BOTH cascades.
    #   16 Auger lines from each scheme, and max|rate difference| / max rate = 0.000e+00 EXACTLY.
    #
    # Branch b: THE EXACT CROSS-CHECK, and the branch to read first.  Dielectronic recombination IS capture
    #   followed by radiative stabilization, so the two schemes must produce the SAME resonances with the SAME
    #   Auger widths -- they differ only in what happens afterwards.  That is an identity, not an
    #   approximation, and it needs no literature value: the same system is run through
    #   Cascade.DielectronicCaptureScheme and Cascade.DielectronicRecombinationScheme and the Auger rates are
    #   compared one by one.
    #
    #   IT IS A REAL CHECK BECAUSE THE TWO PATHS THROUGH THE CODE ARE DIFFERENT.  The capture scheme builds
    #   three groups of configurations and nine Auger steps; the recombination scheme builds a different block
    #   set and three Auger steps plus nine radiative ones.  Agreement to the last bit therefore says the
    #   resonance identification, the block generation and the continuum orbitals all coincide -- which is a
    #   good deal more than a single scheme checked against itself.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fm-crosscheck.sum")
    #
    wa = Cascade.Computation(Cascade.Computation(); name="KLL capture", nuclearModel=Nuclear.Model(6.),
                             grid=grid, approach=Cascade.AverageSCA(), scheme=captureScheme(),
                             initialConfigs=[Configuration("1s^2")] )
    wc = perform(wa; output=true, outputToFile=false)
    drScheme = Cascade.DielectronicRecombinationScheme([E1], false, Shell("2p"), 500.0, 0., 0., 1, [Shell("1s")],
                                                       [Shell("2s"), Shell("2p")], [Shell("2s"), Shell("2p")],
                                                       [Shell("1s"), Shell("2s"), Shell("2p")])
    wd = Cascade.Computation(Cascade.Computation(); name="KLL DR", nuclearModel=Nuclear.Model(6.),
                             grid=grid, approach=Cascade.AverageSCA(), scheme=drScheme,
                             initialConfigs=[Configuration("1s^2")] )
    we = perform(wd; output=true, outputToFile=false)
    #
    linesC = wc["dielectronic-capture lines:"]
    linesD = [d for d in we["cascade data:"] if eltype(d.lines) == AutoIonization.Line][1].lines
    println("\n  Auger lines from the CAPTURE scheme:        ", length(linesC))
    println("  Auger lines from the RECOMBINATION scheme:  ", length(linesD))
    if  length(linesC) == length(linesD)
        rc = sort([l.totalRate for l in linesC]);   rd = sort([l.totalRate for l in linesD])
        println("  " * @sprintf("max |rate difference| / max rate = %.3e     (must be 0)",
                                maximum(abs.(rc .- rd)) / maximum(rd)))
    else
        println("  DIFFERENT COUNTS -- the two schemes disagree about which resonances exist.")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  23-Aug-2026
    #   [PROVENANCE: as branch a.]  ~1.3 s warm.
    #   9 Auger steps are set up and only 3 of them produce a line: steps 1, 4 and 7 give 1, 7 and 8 lines,
    #   totalling the 16 of branch a.  The other SIX produce nothing at all.
    #
    # Branch c: WHY SIX OF THE NINE STEPS YIELD NOTHING, which is the same fact that decides whether a
    #   resonant ionization route exists at all (see example-Fl.jl branch a) and is worth meeting twice.
    #
    #   The capture scheme sets up Auger steps from every resonance to every target block -- to the GROUND
    #   state, which is the capture channel itself, and to the EXCITED states, which is what would turn a
    #   capture into a contribution to electron-impact EXCITATION.  Here the excited channels are all closed:
    #
    #       A RESONANCE CONVERGES TO THE THRESHOLD IT WAS BUILT ON, FROM BELOW.
    #
    #   The 1s 2l 2l' resonances are built on the 1s 2l thresholds, so they lie BELOW them however the capture
    #   shell is chosen, and cannot autoionize into 1s 2s or 1s 2p at all.  Only the channel down to the 1s^2
    #   ground state is open.  That is also exactly why branch b's cross-check comes out identical: with the
    #   excited channels closed, the capture scheme has nothing left that the recombination scheme does not
    #   also compute.
    #
    #   WHAT WOULD OPEN THEM is a resonance built on a HIGHER threshold than the exit channel -- 1s 3l nl'
    #   decaying into 1s 2l, say.  That costs a larger configuration set and is left to a later branch; the
    #   point here is that the closure is structural and not an artefact of the shells chosen.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fm-thresholds.sum")
    #
    wa = Cascade.Computation(Cascade.Computation(); name="KLL capture, step accounting",
                             nuclearModel=Nuclear.Model(6.), grid=grid, approach=Cascade.AverageSCA(),
                             scheme=captureScheme(), initialConfigs=[Configuration("1s^2")] )
    wb = perform(wa; output=true, outputToFile=false)
    linesC = wb["dielectronic-capture lines:"]
    #
    println("\n  Every capture line ends on the SAME target level if the excited channels are closed:")
    finals = unique([ (l.finalLevel.index, l.finalLevel.J, l.finalLevel.parity) for l in linesC ])
    for f in finals
        n = count(l -> l.finalLevel.index == f[1], linesC)
        println("     final level $(f[1])  $(LevelSymmetry(f[2], f[3]))   carries $n of the $(length(linesC)) lines")
    end
    println("\n     One distinct final level means every open channel goes to the ground state of the target,")
    println("     i.e. the six excited-channel steps contributed nothing.  Read the step list printed above:")
    println("     nine Auger steps were set up and only three reported lines.")
    #
    setDefaults("print summary: close", "")
    #
end
