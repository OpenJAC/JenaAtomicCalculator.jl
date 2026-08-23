
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
                                scheme=Cascade.PhotoIonizationScheme([E1], [80.0], Float64[], [Shell("2p")], Shell[],
                                                                     LevelSelection(), [0,1,2], 0., 0.),
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
`TestFrames.testModule_Cascade_Simulation(; short::Bool=true)`  ... tests on module Cascade.
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
