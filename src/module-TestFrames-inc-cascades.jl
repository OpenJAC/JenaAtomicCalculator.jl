
## Functions in this file cover: cascade computations and simulations.
## Alphabetical order within this file.


"""
`TestFrames.testModule_Cascade_PhotonExcitation(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_PhotonExcitation(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("unit: energy", "eV")                       ## the scheme's photon window is read in this unit
    Defaults.setDefaults("print summary: open", "test-Cascade-PhotonExcitation-new.sum")
    printstyled("\n\nTest the module  Cascade for the PhotoExcitationScheme ... \n", color=:cyan)
    ### Make the tests
    ## This test used to have an EMPTY body and returned success = true unconditionally, which is why it did
    ## not notice that every photo-excitation cascade was broken for two days.  It now runs the reference case
    ## of examples/example-Fc.jl: Ne^+ excited from both L subshells into 3s and 3p, giving 2 steps and 27 E1
    ## lines in ~6 s.  The comparison it once carried was disabled because the .jld filename carries a run
    ## date; outputToFile=false removes that filename from the output, so a real comparison is possible again.
    name = "Ne^+ 2s,2p -> 3s,3p photo-excitation"
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, rbox = 11.0)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                approach=Cascade.AverageSCA(),
                                scheme=Cascade.PhotoExcitationScheme([E1], 1.0, 200.0, 1, [Shell("2s"), Shell("2p")],
                                                                     [Shell("3s"), Shell("3p")], LevelSelection(), [0,1], 0., 0.),
                                initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)     ## printing the computation is part of the test: a broken Base.show for a scheme is a real defect
    wb = perform(wa; output=true, outputToFile=false)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-PhotonExcitation-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-PhotonExcitation-new.sum"),
                                "Photoexcitation resonance strength as derived", 15)
    testPrint("testModule_Cascade-PhotonExcitation()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_PhotoAbsorption(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_PhotoAbsorption(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("unit: energy", "eV");   Defaults.setDefaults("unit: cross section", "Mbarn")
    Defaults.setDefaults("print summary: open", "test-Cascade-PhotoAbsorption-new.sum")
    printstyled("\n\nTest the module  Cascade for the PhotoAbsorptionScheme ... \n", color=:cyan)
    ### Make the tests
    ## Beryllium-like carbon near its 1s -> 2p resonance; the same physics as examples/example-Fe.jl but on a
    ## four-electron ion, which brings the computation down to ~8 s.  BOTH steps are exercised: the cascade
    ## computation and the subsequent simulation, the latter because the resonance scale factor of
    ## PhotoExcitation.estimateCrossSection lives there and was once wrong by 740x.  The computed results are
    ## handed to the simulation in memory, wrapped exactly as JLD2.load() would return them, so that no file
    ## with a run-date in its name enters the comparison.
    grid = Radial.Grid(Radial.Grid(false); rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    wa   = Cascade.Computation(Cascade.Computation(); name="Photoabsorption of Be-like C", grid=grid,
                                nuclearModel=Nuclear.Model(6.), approach=Cascade.AverageSCA(),
                                scheme=Cascade.PhotoAbsorptionScheme([E1], [300.0], Float64[],
                                                                     [Shell("1s"), Shell("2s")], [Shell("2p")],
                                                                     LevelSelection(), [0,1], true, true, 0., 0.),
                                initialConfigs=[Configuration("1s^2 2s^2")] )
    println(wa)     ## printing the computation is part of the test, see testModule_Cascade_PhotonIonization
    wb = perform(wa; output=true, outputToFile=false)
    #
    property = Cascade.PhotoAbsorptionSpectrum(true, true, 0.2, 1.0, [en for en = 290.0:0.2:296.0],
                                               Shell[], [(1, 1.0)], Configuration[])
    wc = Cascade.Simulation(Cascade.Simulation(); name="Photoabsorption of Be-like C: Simulation",
                            computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                            property=property, settings=Cascade.SimulationSettings(false, false, 0.) )
    wd = perform(wc; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-PhotoAbsorption-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-PhotoAbsorption-new.sum"),
                                "Absorption cross sections are determined", 12)
    testPrint("testModule_Cascade-PhotoAbsorption()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_PhotonIonization(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_PhotonIonization(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("print summary: open", "test-Cascade-PhotonIonization-new.sum")
    printstyled("\n\nTest the module  Cascade for the PhotonIonizationScheme ... \n", color=:cyan)
    ### Make the tests
    ## The former case was named "Photoionization of Si-" but computed neon, and it gave BOTH photonEnergies
    ## = [0.5] and electronEnergies = [4.0] -- the very combination that PhotoIonizationScheme's own guard
    ## forbids, and at a photon energy far below the ~40 eV needed to ionize Ne^+ at all.  It also never
    ## printed the computation, which is why a broken Base.show for this scheme survived a green suite.
    ## Replaced by the reference case of examples/example-Fd.jl: Ne^+ ionized out of 2p at 80 eV, ~7 s.
    Defaults.setDefaults("unit: energy", "eV")
    name = "Ne^+ 2p photo-ionization at 80 eV"
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, rbox = 11.0)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                                ## LevelSelection(true, indices=[1]) EXPLICITLY, since 24-Aug-2026.  This test asked for
                                ## LevelSelection(), i.e. an INACTIVE selection meaning every initial level -- and its
                                ## approved file records five lines, the ground level alone.  It got that only because
                                ## Cascade.computeSteps ignored the field and hard-wired index 1; when the field was
                                ## honoured again the test computed both levels of the 2p^5 hole, ten lines, and failed.
                                ## The reference is right and the request was wrong: the ground-level computation is
                                ## what has always been compared here, so it is now asked for rather than obtained by
                                ## accident.  No approved file is re-approved.
                                scheme=Cascade.PhotoIonizationScheme([E1], [80.0], Float64[], [Shell("2p")], Shell[],
                                                                     LevelSelection(true, indices=[1]), [0,1,2], 0., 0.),
                                initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)     ## printing the computation is part of the test, see above
    wb = perform(wa; output=true, outputToFile=false)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data; anchored on the SUMMED cross sections, so that the total table
    # is exercised on the cascade path as well as in testModule_PhotoIonization.
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-PhotonIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-PhotonIonization-new.sum"),
                                "Total photoionization cross sections, summed", 14)
    testPrint("testModule_Cascade-PhotonIonization()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_DielectronicCapture(; short::Bool=true)`  ... tests the Cascade module for the
    DielectronicCaptureScheme.
"""
function testModule_Cascade_DielectronicCapture(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("unit: energy", "eV")
    Defaults.setDefaults("print summary: open", "test-Cascade-DielectronicCapture-new.sum")
    printstyled("\n\nTest the module  Cascade for the DielectronicCaptureScheme ... \n", color=:cyan)
    ### Make the tests
    ## The KLL group of helium-like carbon, as in examples/example-Fm.jl.  ~1.3 s warm, the cheapest cascade
    ## test in the file.
    ##
    ## THIS ONE CARRIES AN EXACT INTERNAL CHECK, which is unusual here and is the reason to prefer it over a
    ## plain smoke run.  Dielectronic recombination IS capture followed by radiative stabilization, so the two
    ## schemes must produce the SAME resonances with the SAME Auger widths; they differ only in what happens
    ## afterwards.  Both are run and the rates compared one by one, and the two reach that point along
    ## genuinely different paths through the code -- the capture scheme builds three configuration groups and
    ## nine Auger steps, the recombination scheme a different block set with three Auger and nine radiative
    ## steps.  Agreement to the last bit therefore constrains the resonance identification, the block
    ## generation and the continuum orbitals at once.
    ##
    ## The approved comparison remains a photograph of this code's own output and establishes no physics; the
    ## capture-against-recombination identity above is the part that could actually fail for a real reason.
    grid   = Radial.Grid(Radial.Grid(false); rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 8.0)
    scheme = Cascade.DielectronicCaptureScheme(500.0, 0., 1, [Shell("1s")], [Shell("2s"), Shell("2p")],
                                                [Shell("2s"), Shell("2p")])
    wa     = Cascade.Computation(Cascade.Computation(); name="KLL dielectronic capture of He-like C", grid=grid,
                                  nuclearModel=Nuclear.Model(6.), approach=Cascade.AverageSCA(), scheme=scheme,
                                  initialConfigs=[Configuration("1s^2")] )
    println(wa)     ## printing the computation is part of the test, see testModule_Cascade_PhotonIonization
    wb = perform(wa; output=true, outputToFile=false)
    #
    drScheme = Cascade.DielectronicRecombinationScheme([E1], false, Shell("2p"), 500.0, 0., 0., 1, [Shell("1s")],
                                                       [Shell("2s"), Shell("2p")], [Shell("2s"), Shell("2p")],
                                                       [Shell("1s"), Shell("2s"), Shell("2p")])
    wc = Cascade.Computation(Cascade.Computation(); name="KLL dielectronic recombination of He-like C", grid=grid,
                              nuclearModel=Nuclear.Model(6.), approach=Cascade.AverageSCA(), scheme=drScheme,
                              initialConfigs=[Configuration("1s^2")] )
    wd = perform(wc; output=true, outputToFile=false)
    #
    linesC = wb["dielectronic-capture lines:"]
    linesD = [d for d in wd["cascade data:"] if eltype(d.lines) == AutoIonization.Line][1].lines
    ## THE LAYOUT HERE IS DICTATED BY testCompareFiles, which compares lines iold+2 ... iold+noLines: it SKIPS
    ## the first line after the anchor and needs noLines lines to exist beyond it.  An anchor that carries the
    ## result itself therefore has nothing to index and walks off the end of the file, which is what it did.
    ## So the anchor is a bare header, a rule follows (skipped anyway, being under five characters of content),
    ## and the four comparable numbers come after it.
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nFinal = length(unique([l.finalLevel.index  for l in linesC]))
    if  length(linesC) == length(linesD)  &&  length(linesC) > 0
        rc = sort([l.totalRate for l in linesC]);   rd = sort([l.totalRate for l in linesD])
        sd = "    max |rate difference| / max rate:  $(maximum(abs.(rc .- rd)) / maximum(rd))"
    else
        sd = "    DIFFERENT COUNTS -- the two schemes disagree about which resonances exist."
    end
    for  line  in  [ "\n  Capture against recombination",
                     "  ----------------------------------------------------------",
                     "    Auger lines from the capture scheme:        $(length(linesC))",
                     "    Auger lines from the recombination scheme:  $(length(linesD))",
                     sd,
                     "    distinct final levels of the capture lines: $nFinal",
                     "    summed capture rate [a.u.]:  $(sum([l.totalRate  for l in linesC]))" ]
        println(line);   if  printSummary   println(iostream, line)   end
    end
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-DielectronicCapture-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-DielectronicCapture-new.sum"),
                                "Capture against recombination", 6)
    testPrint("testModule_Cascade_DielectronicCapture()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_ResonantIonization(; short::Bool=true)`  ... tests the Cascade module for the RESONANT
    channels of the ElectronIonizationScheme, comparing the computed data with `test-Cascade-ResonantIonization-approved.sum`.
    A success::Bool is returned.
"""
function testModule_Cascade_ResonantIonization(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("unit: energy", "eV")
    Defaults.setDefaults("print summary: open", "test-Cascade-ResonantIonization-new.sum")
    printstyled("\n\nTest the module  Cascade for the RESONANT channels of the ElectronIonizationScheme ... \n", color=:cyan)
    ### Make the tests
    ## Lithium-like carbon: an incident electron is CAPTURED into a doubly-excited resonance, which then sheds two
    ## electrons and leaves the ion one charge state higher.  A cut-down version of examples/example-Fl.jl -- one
    ## capture shell instead of two -- which brings it to ~10 s, the same order as the photo-absorption test.
    ##
    ## THIS IS A SANITY TEST, NOT A TEST OF THE PHYSICS, and the distinction matters here more than usual.  What it
    ## pins down is that the scheme still runs end to end, that the resonances are still FOUND and still lie at the
    ## same energies, and that the strengths and branchings do not move.  It does NOT establish that the numbers are
    ## right: nothing in the suite compares them against an independent source.
    ##
    ## THE SHELL LISTS ARE NOT ARBITRARY.  A resonance converges to the threshold it was built on FROM BELOW, so it
    ## can never autoionize into that same threshold; capture into 2s/2p while exciting 1s -> 2p closes the
    ## sequential route entirely and every strength comes out zero.  The excitation list therefore spans 2p AND 3s
    ## while the capture goes into 3s, so that the resonances sit on the 1s2s3s threshold and can reach 1s2s2p,
    ## which still carries the K hole and autoionizes again.  A test on the closed configuration would pass just as
    ## happily on a column of zeros, which is exactly the kind of check this file should not contain.
    ##
    ## Both steps are exercised, the computation and the simulation, and the results are handed over IN MEMORY,
    ## wrapped as JLD2.load() would return them, so that no file with a run-date in its name enters the comparison.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 10.0)
    scheme = Cascade.ElectronIonizationScheme(Float64[], [Shell("1s")], [Shell("2p"), Shell("3s")], collect(0:3), 1, 0.,
                                              Basics.AbstractProcess[ResonantImpactIonization.SequentialAuger(),
                                                                     ResonantImpactIonization.SimultaneousAuger()],
                                              [Shell("3s")], 0.)
    wa   = Cascade.Computation(Cascade.Computation(); name="Resonant ionization of Li-like C", grid=grid,
                                nuclearModel=Nuclear.Model(6.), approach=Cascade.AverageSCA(), scheme=scheme,
                                initialConfigs=[Configuration("1s^2 2s")] )
    println(wa)     ## printing the computation is part of the test, see testModule_Cascade_PhotonIonization
    wb = perform(wa; output=true, outputToFile=false)
    #
    ## dblAugerProbability is set to a definite 0.05 rather than left at its 0. default, so that the simultaneous
    ## column carries numbers and a regression in it would be visible; the value itself is arbitrary and the
    ## property's docstring says why the cascade will not choose one for you.
    property = Cascade.ResonantIonizationStrengths(1, 0., 0.05)
    wc = Cascade.Simulation(Cascade.Simulation(); name="Resonant ionization of Li-like C: Simulation",
                            computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                            property=property, settings=Cascade.SimulationSettings(false, false, 0.) )
    wd = perform(wc; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-ResonantIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-ResonantIonization-new.sum"),
                                "Resonance strengths of the resonant electron-capture channels", 22)
    testPrint("testModule_Cascade_ResonantIonization()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_EiiRateCoefficients(; short::Bool=true)`
    ... tests the Cascade.EiiRateCoefficients simulation property, i.e. the fold of the resonant ionization strengths
        with a Maxwellian.  A success::Bool is returned.
"""
function testModule_Cascade_EiiRateCoefficients(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")
    Defaults.setDefaults("method: normalization, pure sine")
    Defaults.setDefaults("unit: energy", "eV")
    Defaults.setDefaults("print summary: open", "test-Cascade-EiiRateCoefficients-new.sum")
    printstyled("\n\nTest the module  Cascade for the EII rate coefficients of the ElectronIonizationScheme ... \n", color=:cyan)
    ### Make the tests
    ## The same cut-down Li-like carbon case as testModule_Cascade_ResonantIonization -- see the notes there on why
    ## the shell lists are what they are -- but carried one step further: the energy-integrated strengths are folded
    ## with a Maxwellian to give the plasma rate coefficient alpha^EII (T), which is the quantity an ionization
    ## balance consumes and the one that can be compared with other codes.
    ##
    ## WHAT THIS ADDS OVER THE STRENGTH TEST, stated narrowly because it is easy to overclaim, and this comment has
    ## already been wrong twice.  It pins the Maxwellian FOLD -- the exponent, the sign of the exponential, the
    ## Kelvin-to-Hartree conversion -- and it does so through the APPROVED NUMBERS, which change wholesale if any of
    ## those is wrong.  That is where the regression power lives.
    ##
    ## It does NOT check that the resonances were correctly identified: prediction and curve are built from the same
    ## energies, so a misplaced set moves both together and passes.
    ##
    ## THE SIX TEMPERATURES ARE CHOSEN, NOT ROUNDED OFF.  The printed bracket check -- is the tabulated maximum one
    ## of the two grid points around kT = 2E/3? -- can only discriminate if the grid resolves the peak.  Measured on
    ## this system, E = 306.6 eV and the maximum sits at 204.4 eV: with the four temperatures this test first used,
    ## replacing T^(-3/2) by T^(-1/2) left the argmax at 430.9 eV, still inside the bracket, so the diagnostic
    ## passed on broken arithmetic.  With 2.0e6 K added the argmax moves to 172.3 eV for the correct exponent and to
    ## 861.7 eV for the wrong one, which is outside the bracket and fails.  A diagnostic that cannot fail on the grid
    ## it is printed for is not worth printing.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 10.0)
    scheme = Cascade.ElectronIonizationScheme(Float64[], [Shell("1s")], [Shell("2p"), Shell("3s")], collect(0:3), 1, 0.,
                                              Basics.AbstractProcess[ResonantImpactIonization.SequentialAuger(),
                                                                     ResonantImpactIonization.SimultaneousAuger()],
                                              [Shell("3s")], 0.)
    wa   = Cascade.Computation(Cascade.Computation(); name="Resonant ionization of Li-like C", grid=grid,
                                nuclearModel=Nuclear.Model(6.), approach=Cascade.AverageSCA(), scheme=scheme,
                                initialConfigs=[Configuration("1s^2 2s")] )
    wb = perform(wa; output=true, outputToFile=false)
    #
    property = Cascade.EiiRateCoefficients(1, [3.0e5, 1.0e6, 2.0e6, 5.0e6, 1.0e7, 3.0e7], 0., 0.05)
    wc = Cascade.Simulation(Cascade.Simulation(); name="Resonant ionization of Li-like C: rate coefficients",
                            computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                            property=property, settings=Cascade.SimulationSettings(false, false, 0.) )
    wd = perform(wc; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-EiiRateCoefficients-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-EiiRateCoefficients-new.sum"),
                                "Electron-impact ionization plasma rate coefficients", 10)
    testPrint("testModule_Cascade_EiiRateCoefficients()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_Simulation(; short::Bool=true)`  ... tests the Cascade module for SIMULATIONS, deriving
    the properties from stored stepwise-decay data (`test/approved/test-Cascade-StepwiseDecay-data.jld`) rather than
    recomputing the cascade. A success::Bool is returned.
"""
function testModule_Cascade_Simulation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Cascade-Simulation-new.sum")
    printstyled("\n\nTest the module  Cascade for Simulations ... \n", color=:cyan)
    ### Make the tests
    datafile = joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-StepwiseDecay-data.jld")
    data = [JLD2.load(datafile)]
    name = "Simulation of the neon 1s^-1 3p decay"

    ## The initial occupations moved from Cascade.SimulationSettings onto Cascade.IonDistribution itself;
    ## SimulationSettings now carries only (printTree, printLongTree, initialPhotonEnergy).  The six-argument
    ## form this test used no longer exists, which is one of the three reasons it had stopped running.
    wc   = Cascade.Simulation(Cascade.Simulation(), name=name,
                                property=Cascade.IonDistribution([(1, 2.0), (2, 1.0), (3, 0.5)], Configuration[]),
                                settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    wd = perform(wc; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-Simulation-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-Simulation-new.sum"), "(Final) Ion distribution for", 8)
    testPrint("testModule_Cascade-Simulation()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_StepwiseDecay(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_StepwiseDecay(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("print summary: open", "test-Cascade-StepwiseDecay-new.sum")
    printstyled("\n\nTest the module  Cascade for the StepwiseDecayScheme ... \n", color=:cyan)
    ### Make the tests
    name = "Cascade after neon 1s --> 3p excitation"
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, rbox = 9.5)
    decayShells = [Shell(1,0), Shell(2,0), Shell(2,1), Shell(3,1)]
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                                scheme=Cascade.StepwiseDecayScheme([Auger(), Radiative()], 1, Dict{Int64,Float64}(), 0, decayShells, Shell[], Shell[]),
                                initialConfigs=[Configuration("1s^1 2s^2 2p^6 3p")] )
    println(wa)
    wb = perform(wa; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-StepwiseDecay-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-StepwiseDecay-new.sum"), "Steps that are defined for the curren", 20)
    testPrint("testModule_Cascade-StepwiseDecay()::", success)
    return(success)
end
