
println("Rb) Apply & test the NistCalibration module: replace computed level energies by measured ones and")
println("    show the effect on a real computed rate -- not just a nicer-looking energy table.")

setDefaults("print summary: open", "zzz-NistCalibration.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

if  true
    # Last visit:  18-Oct-2026
    # Last successful:  14-Oct-2026
    #
    # Branch a: Na I (Z=11), 3s ^2S_1/2 (ground) -> 3p ^2P_1/2, ^2P_3/2 -- the sodium D-lines, among the most
    #   precisely measured transitions in existence, and already a system this project's own examples use
    #   (see example-Cb.jl's Na HFS branch, where the field choice alone moves a computed quantity by 26% --
    #   Na's single valence electron outside a polarisable Ne-like core is NOT hydrogenic, and a plain
    #   single-configuration ab initio treatment is known to be noticeably off here).
    #
    # NIST VALUES QUOTED, NOT RE-FETCHED. 3p_1/2 at 16956.170 cm^-1 and 3p_3/2 at 16973.379 cm^-1 above the
    # 3s ground are well-established spectroscopic constants (the D1/D2 line origins), consistent across NIST
    # ASD and essentially every atomic-physics reference -- quoted here from general knowledge, not a literal
    # downloaded snapshot, and the branch says so honestly rather than inventing a specific file date.
    ni = Nuclear.Model(11.)
    wa = Atomic.Computation(Atomic.Computation(), name="Rb-a: Na I, 3s -> 3p", grid=Radial.Grid(true),
                            nuclearModel=ni,
                            configs=[Configuration("[Ne] 3s"), Configuration("[Ne] 3p")],
                            propertySettings=[ Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000.) ] )
    wb    = perform(wa; output=true)
    lines = wb[Einstein.Settings]
    #
    # The RAW (ab initio) A-values, straight from the lines Einstein already computed:
    println("\n>>> RAW (ab initio) A-values:")
    for  ln in lines
        ln.initialLevel.energy <= ln.finalLevel.energy   &&   continue
        deltaE = Defaults.convertUnits("energy: from atomic to Kayser", ln.initialLevel.energy - ln.finalLevel.energy)
        aVal   = Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin)
        @printf("    J_i=%s -> J_f=%s   deltaE = %10.3f cm^-1 (ab initio)   A(Bab) = %.5e /s\n",
                ln.initialLevel.J, ln.finalLevel.J, deltaE, aVal)
    end
    #
    # Build the multiplet NistCalibration needs directly from the SAME levels Einstein already touched --
    # ground (3s, lowest energy among the lines' own levels) and the two 3p levels.
    allLevels = unique(vcat([ln.initialLevel for ln in lines], [ln.finalLevel for ln in lines]))
    ground    = allLevels[ argmin([lv.energy for lv in allLevels]) ]
    p12       = [lv for lv in allLevels if lv.J == AngularJ64(1//2) && lv.index != ground.index][1]
    p32       = [lv for lv in allLevels if lv.J == AngularJ64(3//2)][1]
    mp        = Multiplet("Na I 3s,3p", [ground, p12, p32])

    nistTable = [ NistCalibration.NistEntry(p12.index, 16956.170, "3p_1/2 (D1 line origin)"),
                  NistCalibration.NistEntry(p32.index, 16973.379, "3p_3/2 (D2 line origin)") ]
    cm = NistCalibration.calibrate(mp, nistTable,
                     "NIST ASD, well-established D-line constants (not a dated snapshot -- see this branch's own comment)";
                     groundLevelIndex=ground.index)
    NistCalibration.displayResults(stdout, cm)
    println("\n>>> ", NistCalibration.summary(cm))
    #
    # NOW RECOMPUTE THE A-VALUES FROM THE CALIBRATED MULTIPLET -- same wavefunctions, calibrated energies.
    # This is the point of the whole module: cm.multiplet is an ORDINARY Multiplet, so it goes straight into
    # the SAME Einstein.computeLines JAC would always use, no change to Einstein needed at all.
    calibLines = Einstein.computeLines(cm.multiplet, wa.grid,
                                        Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000.))
    println("\n>>> CALIBRATED A-values (same matrix elements, NIST-calibrated energies):")
    for  ln in calibLines
        ln.initialLevel.energy <= ln.finalLevel.energy   &&   continue
        deltaE = Defaults.convertUnits("energy: from atomic to Kayser", ln.initialLevel.energy - ln.finalLevel.energy)
        aVal   = Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin)
        @printf("    J_i=%s -> J_f=%s   deltaE = %10.3f cm^-1 (NIST)        A(Bab) = %.5e /s\n",
                ln.initialLevel.J, ln.finalLevel.J, deltaE, aVal)
    end
    println("\n>>> Report this table with its own provenance line, since nothing downstream carries it:")
    println(">>> ", NistCalibration.summary(cm))
    #
end
