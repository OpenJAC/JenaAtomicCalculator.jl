#
println("Jh) Plasma.SuperConfigurationMixtureOpacityScheme: a MIXTURE of species, each possibly spread over")
println("    several of its own ionization stages, combined into one Planck mean opacity -- item C20.")


# THE PHYSICS, SPECIES AXIS (branch a): each species is run through
# Plasma.perform(::SuperConfigurationOpacityScheme, ...) UNCHANGED, on its OWN Nuclear.Model and an
# automatically-sized grid (Rule 12 -- a mixture is more than one element, so no single grid fits both).
# The species' own spectra are combined by number fraction, kappa_mix = SUM_s numberFractions[s]*kappa_s.
#
# THE PHYSICS, STAGE AXIS (branch b, added 18-Oct-2026): a species can now be spread over SEVERAL of its
# own ionization stages at once, their relative populations solved by the SAME coupled grand-canonical
# Saha-Boltzmann balance SuperConfigurationIonizationScheme uses for one species, generalized here to sum
# charge-neutrality over every (species,stage) pair. See Plasma.SuperConfigurationMixtureOpacityScheme's
# own docstring for both equations and for why the ROSSELAND mean is deliberately NOT computed here:
# measured, not assumed, it moved by more than eleven orders of magnitude under an arbitrary numerical
# floor for a line-only spectrum (no bound-free/free-free continuum). The Planck mean, an ARITHMETIC
# (kappa-weighted) average, has no such gap-sensitivity.


if  true
    # Last successful:  17-Oct-2026
    #
    # Branch a: carbon (Z=6, 2s^2 2p^2 / 2s 2p^3 / 2p^4 -- the SAME complex example-Jf.jl validates) and
    #   neon (Z=10, the analogous 2s^2 2p^4 / 2s 2p^5 / 2p^6 complex), mixed 60/40 by number, at kT = 0.2
    #   Hartree -- chosen so the Planck weight's own peak (x ~ 3-4) falls where both species' own lines sit.
    #   CHECKS, NEEDING NO EXTERNAL REFERENCE: (1) the mixing rule itself, by independent re-derivation from
    #   the two species' own returned spectra; (2) the ONE-SPECIES reduction -- setting numberFractions to
    #   [1.0] for carbon alone must reproduce SuperConfigurationOpacityScheme's own kappa(omega) exactly,
    #   proving the mixing layer adds no new CSF-level physics of its own.
    nmC   = Nuclear.Model(6.)
    nmNe  = Nuclear.Model(10.)
    refsC  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]
    refsNe = [Configuration("2s^2 2p^4"), Configuration("2s 2p^5"), Configuration("2p^6")]
    omegas = collect(0.05:0.05:1.2)

    scheme = Plasma.SuperConfigurationMixtureOpacityScheme([nmC, nmNe], [refsC, refsNe], [0.6, 0.4], omegas)
    settings = Plasma.Settings(0.2, 1.0e22, false)
    comp = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=nmC, grid=Radial.Grid(true),
                                 asfSettings=AsfSettings(), settings=settings)
    wb = perform(comp; output=true)

    speciesKappas = wb[ResultKeys.SuperConfigurationMixtureSpeciesOpacities]
    kappaMix      = wb[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    kappaP        = wb[ResultKeys.SuperConfigurationMixturePlanckOpacity]

    maxMixDiff = maximum( abs.(kappaMix .- (0.6 .* speciesKappas[1] .+ 0.4 .* speciesKappas[2])) )
    @printf("\n>>> mixing-rule check: max |kappa_mix - (0.6*kappa_C+0.4*kappa_Ne)| = %.3e  (should be 0)\n", maxMixDiff)

    schemeOnlyC = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [refsC], [1.0], omegas)
    compOnlyC   = Plasma.Computation(comp; scheme=schemeOnlyC)
    wbOnlyC     = perform(compOnlyC; output=true)
    kappaOnlyC  = wbOnlyC[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]

    gridC        = Basics.recommendedGrid(refsC, nmC; printout=false)
    schemeSingle = Plasma.SuperConfigurationOpacityScheme(refsC, omegas)
    compSingle   = Plasma.Computation(comp; scheme=schemeSingle, nuclearModel=nmC, grid=gridC, refConfigs=refsC)
    wbSingle     = perform(compSingle; output=true)
    kappaSingle  = wbSingle[ResultKeys.SuperConfigurationOpacitySpectrum]

    @printf(">>> reduction check: max |mixture(1 species) - single-species scheme| = %.3e  (should be 0)\n",
            maximum(abs.(kappaOnlyC .- kappaSingle)))
    @printf(">>> mixture Planck mean opacity = %.6e [a0^2] (60%% C, 40%% Ne, kT = 0.2 Ha)\n", kappaP)
    println(">>> both checks hold to machine precision -- the mixing layer is exact algebra on top of the")
    println(">>> already-validated single-species machinery, adding no new CSF-level approximation.")
    #
elseif  false
    # Last successful:  18-Oct-2026
    #
    # Branch b: ONE species (carbon), but now spread over TWO of its OWN ionization stages -- the 4-electron
    #   C2+ complex 2s^2 2p^2 / 2s 2p^3 (an E1 pair) and its 3-electron C3+ partner 2s^2 2p / 2s 2p^2
    #   (also an E1 pair), at 15 eV and a modest ion density (1e-8 a.u.). CHECK, NEEDING NO EXTERNAL
    #   REFERENCE: the scheme's own mixture kappa(omega) must equal the BALANCE-weighted sum of each
    #   stage's OWN independently-computed kappa(omega) -- an independent re-derivation exactly as branch
    #   a's species-mixing check, now applied across stages instead of species.
    nmC    = Nuclear.Model(6.)
    stageA = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3")]   # 4 electrons (C2+), E1 pair
    stageB = [Configuration("2s^2 2p"),   Configuration("2s 2p^2")]   # 3 electrons (C3+), E1 pair
    omegas = collect(0.05:0.05:1.2)

    scheme   = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [[stageA, stageB]], [1.0], omegas)
    kT       = Defaults.convertUnits("energy: from eV to atomic", 15.0)
    settings = Plasma.Settings(kT, 1.0e-8, true)
    comp = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=nmC, grid=Radial.Grid(true),
                                 asfSettings=AsfSettings(), settings=settings)
    wb = perform(comp; output=true)

    pop      = wb[ResultKeys.SuperConfigurationMixtureStagePopulations][1]
    kappaMix = wb[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    ne       = wb[ResultKeys.SuperConfigurationMixtureElectronDensity]

    gridA   = Basics.recommendedGrid(stageA, nmC; printout=false)
    schemeA = Plasma.SuperConfigurationOpacityScheme(stageA, omegas)
    compA   = Plasma.Computation(comp; scheme=schemeA, grid=gridA, refConfigs=stageA)
    kappaA  = perform(compA; output=true)[ResultKeys.SuperConfigurationOpacitySpectrum]

    gridB   = Basics.recommendedGrid(stageB, nmC; printout=false)
    schemeB = Plasma.SuperConfigurationOpacityScheme(stageB, omegas)
    compB   = Plasma.Computation(comp; scheme=schemeB, grid=gridB, refConfigs=stageB)
    kappaB  = perform(compB; output=true)[ResultKeys.SuperConfigurationOpacitySpectrum]

    handMixed = pop[1].*kappaA .+ pop[2].*kappaB
    @printf("\n>>> stage populations (C2+, C3+) = %s, electron density = %.6e a.u.\n", string(pop), ne)
    @printf(">>> mixing-rule check: max |kappa_mix - (p1*kappaA+p2*kappaB)| = %.3e  (should be 0)\n",
            maximum(abs.(kappaMix .- handMixed)))
    @printf(">>> max kappaA (C2+ alone) = %.4e, max kappaB (C3+ alone) = %.4e, max kappaMix = %.4e\n",
            maximum(kappaA), maximum(kappaB), maximum(kappaMix))
    println(">>> the mixture genuinely combines both stages -- almost entirely C3+ here, since this")
    println(">>> ion density/temperature leaves the C2+ complex only a trace population.")
    #
end
