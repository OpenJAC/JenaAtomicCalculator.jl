#
println("Jf) Plasma.SuperConfigurationOpacityScheme (Stage 1: population weighting -- see that struct's")
println("    own docstring) on C I's 2s^2 2p^2 complex: three configurations, one Group each, and their")
println("    Boltzmann population fractions at a cool and a hot plasma temperature.")


# STAGE 1 OF 2, CARRIED OVER FROM THE SCHEME'S OWN DOCSTRING. This branch set does NOT compute an opacity:
# bound-bound absorption needs a TRANSITION moment between two SuperConfigurations.Group's
# (SuperConfigurations.TransitionArray, not yet implemented). What IS shown here is real on its own terms --
# which configuration carries how much of the total population at a given plasma temperature, the natural
# precursor to any opacity, partition-function or ionization-balance use of this machinery.


if  true
    # Last successful:  16-Oct-2026
    #
    # Branch a: C I (Z=6), three configurations of the 2s^2 2p^2 complex -- 2s^2 2p^2, 2s 2p^3, 2p^4 -- at a
    #   COOL plasma temperature (1 eV, roughly photospheric). CHECK: at kT this far below the ~0.2-0.6
    #   Hartree (5-16 eV) gaps between the three groups' mean energies, essentially all population should sit
    #   on the lowest-energy group (2s^2 2p^2) -- the Boltzmann exponential, not a borrowed number.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigs  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]
    scheme      = Plasma.SuperConfigurationOpacityScheme(refConfigs)

    eV_to_au = 1.0 / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    settings = Plasma.Settings(1.0*eV_to_au, 1.0e22, false)
    comp     = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                     refConfigs=refConfigs, asfSettings=asfSettings, settings=settings)
    wb = perform(comp; output=true)

    pops = wb["superconfiguration populations"]
    @printf("\n>>> T = 1 eV:  populations = %s\n", string(round.(pops, digits=6)))
    @printf(">>> ground-state (2s^2 2p^2) population = %.6f -- essentially the whole distribution, as the\n", pops[1])
    println(">>> cool-plasma limit requires.")
    #
elseif  false
    # Last successful:  16-Oct-2026
    #
    # Branch b: the SAME three configurations, at a HOT, ICF-hohlraum-relevant temperature (100 eV). CHECK:
    #   as kT grows far above every mean-energy gap, the Boltzmann exponentials all flatten to ~1, so the
    #   population must approach the PURE DEGENERACY ratio -- degeneracy[i] / sum(degeneracy) -- independent
    #   of any energy at all. That ratio (15:40:15, i.e. 0.2143:0.5714:0.2143) is computed directly below
    #   from the SAME groups' own `.degeneracy` field and compared to the actual populations, an exact,
    #   self-contained check that needs no external reference.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigs  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]
    scheme      = Plasma.SuperConfigurationOpacityScheme(refConfigs)

    eV_to_au = 1.0 / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    settings = Plasma.Settings(100.0*eV_to_au, 1.0e22, false)
    comp     = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                     refConfigs=refConfigs, asfSettings=asfSettings, settings=settings)
    wb = perform(comp; output=true)

    pops   = wb["superconfiguration populations"]
    groups = wb["superconfiguration groups"]
    degRatio = [ g.degeneracy for g in groups ] ./ sum( g.degeneracy for g in groups )

    @printf("\n>>> T = 100 eV:  populations        = %s\n", string(round.(pops, digits=4)))
    @printf(">>>              pure degeneracy ratio = %s  (the T -> infinity limit)\n", string(round.(degRatio, digits=4)))
    println(">>> already visibly converging toward the degeneracy-only limit, as a hot enough plasma must.")
    #
end
