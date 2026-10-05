
println("Dj) Test of the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoDoubleIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
quasiShells = Basics.generateShellList(4, 4, 1)

if  true
    # Last successful:  unknown ...
    # Photo-double ionization of he-like Be: 1s^2   -->  4l + epsilon kappa: Generation of MeanFieldMultiplet
    name        = "He-like Be"
    refConfigs  = Basics.generateConfigurationsWithAdditionalElectron([Configuration("1s")], quasiShells)
    mfSettings  = MeanFieldSettings()
    #
    wa          = Representation(name, Nuclear.Model(4.01), Radial.Grid(true), refConfigs, MeanFieldMultiplet(mfSettings) )
    wb          = generate(wa, output=true)
    #
elseif true
    # Last successful:  unknown ...
    # Photo-double ionization of he-like Be: 1s^2   -->  4l + epsilon kappa: Computation of cross sections
    # `wb` is produced by an earlier branch of this file (line 21).
    if  !@isdefined(wb)   error("Run the earlier branch of this file first; it defines `wb`.")   end
    nMultiplet     = wb[ResultKeys.MeanFieldMultiplet]
    # The KEYWORD form, so that a new Settings field cannot silently break this call -- which had already happened:
    # the positional version passed `quasiShells` third, where `photonEnergies` now sits, and a field of that name
    # no longer exists.  A second field was added on 30-Sep-2026 (NoIntermediateEnergies) and would have been the
    # eleventh of twelve.
    doubleSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                         multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[400., 500.],
                         NoEnergySharings=2, maxKappa=2, calcDifferentialCs=true, printBefore=true,
                         lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                         eeInteraction=CoulombInteraction(), gMultiplet=nMultiplet )
    
    wc = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(4.01), 
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s^0")],
                            processSettings = doubleSettings )

    wd = perform(wc)
    #
end
#
setDefaults("print summary: close", "")


