
println("Rc) Apply & test the UncertaintyEstimate module: the spread of a real rate across a FEW named,")
println("    physically motivated variants -- a mean-field choice (Type A) and a correlation layer (Type B).")

setDefaults("print summary: open", "zzz-UncertaintyEstimate.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

if  true
    # Last visit:  20-Oct-2026
    # Last successful:  14-Oct-2026
    #
    # Na I (Z=11), 3s ^2S_1/2 -> 3p ^2P_3/2 -- the SAME system and transition example-Rb.jl already used for
    # NistCalibration, reused here rather than introducing a third system for one more demonstration.
    #
    # TYPE A -- MEAN-FIELD CHOICE.  "did you optimize the orbitals for the configuration average (DFS), or
    # would a level-optimized field (AL) answer differently?"  A question every atomic physicist recognises
    # without needing any numerical explanation, and one this project has already measured moving a computed
    # quantity by 26% for Na specifically (example-Cb.jl's HFS branch).
    function  aValue(scField)
        ni = Nuclear.Model(11.)
        wa = Atomic.Computation(Atomic.Computation(), name="Rc: Na I, field=$scField", grid=Radial.Grid(true),
                                nuclearModel=ni,
                                configs=[Configuration("[Ne] 3s"), Configuration("[Ne] 3p")],
                                asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=scField),
                                propertySettings=[ Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000.) ] )
        wb    = perform(wa; output=true)
        lines = wb[Einstein.Settings]
        ln    = [l for l in lines if l.initialLevel.J == AngularJ64(3//2) && l.finalLevel.J == AngularJ64(1//2) &&
                                      l.initialLevel.energy > l.finalLevel.energy][1]
        return( Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin) )
    end

    variantsA = [ UncertaintyEstimate.Variant("DFS", aValue(Basics.DFSField())),
                  UncertaintyEstimate.Variant("AL",  aValue(Basics.ALField())) ]
    estA = UncertaintyEstimate.estimate("A(3p_3/2 -> 3s), Na I -- mean-field choice", variantsA)
    UncertaintyEstimate.displayResults(stdout, estA)

    #
    # TYPE B -- CORRELATION EXTENT.  "does adding one more correlating configuration move the answer?" --
    # the same question Rule 21's layer-increment comparisons already ask by hand, formalised here as a
    # reusable spread rather than a one-off, bespoke comparison. DFS throughout, so only the correlation
    # space (not the field) differs between the two variants.
    function  aValueCorr(configs)
        ni = Nuclear.Model(11.)
        wa = Atomic.Computation(Atomic.Computation(), name="Rc: Na I, corr", grid=Radial.Grid(true),
                                nuclearModel=ni, configs=configs,
                                propertySettings=[ Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000.) ] )
        wb    = perform(wa; output=true)
        lines = wb[Einstein.Settings]
        ln    = [l for l in lines if l.initialLevel.J == AngularJ64(3//2) && l.finalLevel.J == AngularJ64(1//2) &&
                                      l.initialLevel.energy > l.finalLevel.energy][1]
        return( Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin) )
    end

    variantsB = [ UncertaintyEstimate.Variant("ref-only",   aValueCorr([Configuration("[Ne] 3s"), Configuration("[Ne] 3p")])),
                  UncertaintyEstimate.Variant("ref+1layer", aValueCorr([Configuration("[Ne] 3s"), Configuration("[Ne] 3p"), Configuration("[Ne] 3d")])) ]
    estB = UncertaintyEstimate.estimate("A(3p_3/2 -> 3s), Na I -- correlation extent", variantsB)
    UncertaintyEstimate.displayResults(stdout, estB)
    #
end
