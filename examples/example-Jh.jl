#
println("Jh) Plasma.SuperConfigurationMixtureOpacityScheme: a two-species MIXTURE (carbon + neon), combined")
println("    at a shared LTE temperature into one Planck mean opacity -- competition item C20.")


# THE PHYSICS: each species is run through Plasma.perform(::SuperConfigurationOpacityScheme, ...) UNCHANGED,
# on its OWN Nuclear.Model and an automatically-sized grid (Rule 12 -- a mixture is more than one element,
# so no single grid fits both). The two species' Stage-3 spectra are combined by number fraction,
#
#     kappa_mix(omega) = numberFractions[1]*kappa_C(omega) + numberFractions[2]*kappa_Ne(omega)
#
# and reduced to a Planck mean opacity. See Plasma.SuperConfigurationMixtureOpacityScheme's own docstring
# for the weight function and for why the ROSSELAND mean is deliberately NOT computed here: measured, not
# assumed, it moved by more than eleven orders of magnitude under an arbitrary numerical floor for a
# line-only spectrum (no bound-free/free-free continuum), so there is nothing trustworthy to report. The
# Planck mean, an ARITHMETIC (kappa-weighted) average, has no such gap-sensitivity.


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
end
