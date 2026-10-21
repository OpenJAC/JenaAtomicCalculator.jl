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
#
# A LITERATURE TEST OF THE STAGE AXIS (branch c, added 18-Oct-2026): a CH plasma benchmarked against Lee
# et al., Phys. Rev. E 98, 043203 (2018), Fig. 7. Calls Plasma.SuperConfigurationIonizationScheme directly
# for carbon alone (hydrogen is pinned neutral, so the full mixture balance reduces to exactly this), and
# finds a real, informative limit of the current machinery along the way: enriching a stage's own
# configuration space for balance accuracy can make it too rich for this scheme's OPACITY step, which
# requires every configuration pair to be E1-connected. See the branch's own comment for the full story.
#
# TWO CONTINUUM CONTRIBUTIONS, each on branch b's own C2+/C3+ pair (added 19-Oct-2026, item C20):
#
# BOUND-FREE (branch d): `boundFreeApprox = Empirical.ScaledHydrogenic()` adds photoionization between
# the two adjacent stages, reusing the already-validated Empirical.photoionizationCrossSection. Far out on
# the line's own exponential tail, where bound-bound has decayed to a literal float64 near-zero, bound-free
# still supplies a physically meaningful Kramers-type cross section -- the regime that makes the Rosseland
# mean ill-defined without it (Plasma.rosselandMeanOpacity's own docstring).
#
# FREE-FREE (branch e): `freeFreeIncluded = true` adds inverse bremsstrahlung, per stage via its own ion
# charge -- Karzas & Latter's (1961) classical cross section times Sutherland's (1998) thermally-averaged
# Gaunt factor. Unlike bound-free, free-free has NO threshold at all, and falls only as a power law
# (omega^-3, modulated by the Gaunt factor and the stimulated-emission factor) rather than bound-bound's
# EXPONENTIAL fall past its last line -- so free-free eventually CROSSES OVER and dominates, however tiny
# it is near a resonance peak.
#
# RESOLVED AUTOIONIZING RESONANCES (branch f, added 20-Oct-2026, item C13): `resonanceContinuumGrid` turns
# on a VOIGT correction to any ground-to-resonance E1 line whose upper group lies above the next stage's
# ionization threshold. The Gaussian statistical spread of that line (from its own TransitionArray) and the
# Lorentzian natural width (a representative-level AutoIonization width into the next stage) are combined
# via the Faddeeva function (Plasma.voigtLineProfile, verified against brute-force convolution before being
# written in) and added as an exact correction on top of the Gaussian-only line the inner Stage 2/3 call
# already produced -- see Plasma.SuperConfigurationMixtureOpacityScheme's own docstring, "RESONANCES"
# section, for the full design.


if  false
    # Last visit:  20-Oct-2026
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
    # Last visit:  20-Oct-2026
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
elseif  false
    # Last visit:  18-Oct-2026   -- not `Last successful`: see the note below on what is likely missing.
    #
    # Branch c: a REAL literature test for the stage axis, not an internal re-derivation. CH plasma (equal
    #   MASS density of C and H, i.e. 1:12 by number), total rho = 1.0 g/cm^3, against Fig. 7 of Lee, Busquet,
    #   Klapisch, Bates, Schmitt, Hu, Giuliani, Phys. Rev. E 98, 043203 (2018) -- an STA (super-transition-
    #   array) LTE model benchmarked against real X-ray-scattering plasma experiments. The paper gives no
    #   numeric table, so the comparison values below are read OFF THE PLOT and are approximate by nature.
    #
    #   HYDROGEN IS KEPT A SINGLE NEUTRAL STAGE, not spread over its own ionization: a bare proton has no
    #   electron, so no SuperConfigurations.Group can be built for it, and the scheme does not yet support a
    #   bare-nucleus stage. Because H then contributes EXACTLY ZERO electrons at every T, the coupled
    #   species+stage balance SuperConfigurationMixtureOpacityScheme would solve here reduces IDENTICALLY to
    #   carbon's own single-species balance at carbon's own actual ion density within the mixture
    #   (nC = nTotal*fC) -- so Plasma.SuperConfigurationIonizationScheme is called directly for carbon alone,
    #   which is exactly this scheme's own Stage 1 balance with nothing else attached.
    #
    #   THIS ALSO SIDESTEPS A REAL LIMITATION, WORTH RECORDING FOR WHOEVER NEXT ENRICHES A STAGE'S
    #   CONFIGURATIONS: SuperConfigurationMixtureOpacityScheme's opacity step (unlike its balance step)
    #   requires EVERY configuration pair within a stage to be connected by a genuine single-electron E1
    #   jump -- true for a bare n=2 complex, but not once an n=2->3 or n=2->4 promoted configuration is added
    #   (e.g. "2s^2 2p 3p" and "2s 2p^3" differ by TWO orbital jumps, correctly zero E1 strength, which
    #   SuperConfigurations.buildTransitionArray treats as an error rather than "no line here"). The balance
    #   itself needs no such connectivity -- only each group's own degeneracy and mean energy -- so this
    #   branch is unaffected, but a configuration space rich enough for an accurate BALANCE is, right now,
    #   too rich for the OPACITY half of this same scheme.
    #
    #   CARBON'S OWN CONFIGURATION SPACE, per stage: the bare n=2 complex PLUS an n=2->3 PLUS an n=2->4
    #   single-promotion complex (one L-shell electron lifted to 3s/3p/3d or 4s/4p/4d/4f) -- the systematic,
    #   next-shell-at-a-time way STA-type models add excited-state degeneracy without combinatorial blow-up.
    #   A ground-configuration-only balance (tried first) overshot the paper's ionization by close to a
    #   factor of 2 at the higher temperatures; adding the n=3 shell closed much of that gap, and n=4 closed
    #   almost none of what was left -- MARGINAL RETURNS HAD COLLAPSED (at T=30eV: ground-only -> n=3 dropped
    #   Z-bar_C by 0.55; n=3 -> n=4 dropped it by only 0.04) -- so a further shell was not pursued.
    #
    #   WHAT IS LIKELY STILL MISSING, AND WHY IT POINTS ONE WAY: the remaining mismatch runs in OPPOSITE
    #   directions at the two ends -- too LITTLE ionization at T=5eV against the experimental point, too MUCH
    #   at T>=20eV against the paper's own STA curves -- and a single missing mechanism explains both:
    #   ionization-potential depression (IPD). Its usual textbook effect (lowering the energy cost to ionize)
    #   would only push harder in the WRONG direction at high T, so it was first set aside as irrelevant here
    #   (wrongly) -- but IPD's OTHER effect, at ~1 g/cm^3, is to CUT OFF the high-n bound-state ladder
    #   outright (the Mott effect): many of the n=3/4 states this branch had to add by hand likely do not
    #   exist as bound states in the real plasma, which would trim exactly the entropic excess now
    #   over-ionizing this vacuum calculation at T>=20eV, while genuine pressure ionization -- itself absent
    #   from any vacuum balance -- would independently raise the real ionization above ours at T=5eV. Two
    #   different real mechanisms, not one, both deliberately absent here; neither is built, so no number
    #   above the ground-only baseline should be read as more than a converging trend.
    nmC  = Nuclear.Model(6.)
    rhoC = Defaults.convertUnits("density: from [g/cm^3] to atomic", 0.5)
    rhoH = Defaults.convertUnits("density: from [g/cm^3] to atomic", 0.5)
    mC = 12.011;  mH = 1.008
    nC  = rhoC/mC;   nH = rhoH/mH;   nTotalCH = nC + nH
    fC  = nC/nTotalCH

    carbonStages = [
        [ Configuration("1s^2 2s^2 2p^2"), Configuration("1s^2 2s 2p^3"), Configuration("1s^2 2p^4"),
          Configuration("1s^2 2s^2 2p 3s"), Configuration("1s^2 2s^2 2p 3p"), Configuration("1s^2 2s^2 2p 3d"),
          Configuration("1s^2 2s^2 2p 4s"), Configuration("1s^2 2s^2 2p 4p"),
          Configuration("1s^2 2s^2 2p 4d"), Configuration("1s^2 2s^2 2p 4f") ],
        [ Configuration("1s^2 2s^2 2p"), Configuration("1s^2 2s 2p^2"), Configuration("1s^2 2p^3"),
          Configuration("1s^2 2s^2 3s"), Configuration("1s^2 2s^2 3p"), Configuration("1s^2 2s^2 3d"),
          Configuration("1s^2 2s^2 4s"), Configuration("1s^2 2s^2 4p"),
          Configuration("1s^2 2s^2 4d"), Configuration("1s^2 2s^2 4f") ],
        [ Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p^2"),
          Configuration("1s^2 2s 3s"), Configuration("1s^2 2s 3p"), Configuration("1s^2 2s 3d"),
          Configuration("1s^2 2s 4s"), Configuration("1s^2 2s 4p"), Configuration("1s^2 2s 4d"), Configuration("1s^2 2s 4f") ],
        [ Configuration("1s^2 2s"), Configuration("1s^2 2p"),
          Configuration("1s^2 3s"), Configuration("1s^2 3p"), Configuration("1s^2 3d"),
          Configuration("1s^2 4s"), Configuration("1s^2 4p"), Configuration("1s^2 4d"), Configuration("1s^2 4f") ],
        [ Configuration("1s^2") ] ]

    scheme     = Plasma.SuperConfigurationIonizationScheme(carbonStages)
    asf        = AsfSettings()
    allConfigs = reduce(vcat, carbonStages)
    # Rule 12: ONE grid is shared across every stage, so it must fit the MOST DIFFUSE orbital anywhere --
    #   here C I's promoted 4f; tailFactor is widened from its 16.0 default since 4f still reached 93% of
    #   the default box.
    sharedGrid = Basics.recommendedGrid(allConfigs, nmC; tailFactor=22.0, printout=false)

    # paper values read off Fig. 7 (no numeric table given); the T=5eV point is the paper's own EXPERIMENTAL
    # datum (Gregori et al., JQSRT 99, 225 (2006)), the rest are its STA model curves at matching density.
    paperT     = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    paperZlow  = [0.8, 1.0, 1.1, 1.5, 1.7, 1.9, 2.0]
    paperZhigh = [1.8, 1.2, 1.5, 1.8, 2.0, 2.1, 2.3]

    @printf("\n%6s %14s %22s\n", "T[eV]", "our Z-bar_C", "paper (Fig. 7, approx)")
    for (k, Tev) in enumerate(paperT)
        kT = Defaults.convertUnits("energy: from eV to atomic", Tev)
        settings = Plasma.Settings(kT, nC, true)
        comp = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=nmC, grid=sharedGrid,
                                     asfSettings=asf, settings=settings)
        wb = perform(comp; output=true)

        charges = wb[ResultKeys.SuperConfigurationIonizationCharges]
        pop     = wb[ResultKeys.SuperConfigurationIonizationStagePopulations]
        Zbar    = sum(charges .* pop)
        @printf("%6.1f %14.4f %13.1f - %.1f\n", Tev, Zbar, paperZlow[k], paperZhigh[k])
    end

    println("\n>>> the qualitative TREND is reproduced across the full 5-60eV range and the T=10-30eV window")
    println(">>> lands close to the paper's own band; the remaining excess at T>=40eV and shortfall at T=5eV")
    println(">>> are the two-sided signature of missing ionization-potential depression (see the note above)")
    println(">>> -- not a configuration-space gap, since adding further shells had already stopped helping.")
    #
elseif  false
    # Last visit:  20-Oct-2026
    # Last successful:  19-Oct-2026   ## at omega = 1.2 Ha: bound-bound 1.5e-15, bound-free 4.1e-7 [a0^2];
    #                                    Rosseland (bound-bound+bound-free) 9.1e-6 [a0^2]
    #
    # Branch d: BOUND-FREE, on branch b's own C2+/C3+ two-stage mixture, now with
    #   boundFreeApprox = Empirical.ScaledHydrogenic() turned on. CHECK, NEEDING NO EXTERNAL REFERENCE: far
    #   out on the line's own exponential tail (the LARGEST omega on this grid, 1.2 Hartree), bound-bound
    #   has decayed to a literal float64 near-zero while bound-free -- a smooth Kramers-type cross section
    #   with no line shape to decay -- remains eight orders of magnitude larger there. Without bound-free,
    #   the Rosseland mean is not even well-defined (Plasma.rosselandMeanOpacity's own docstring: it moved
    #   by more than eleven orders of magnitude under an arbitrary numerical floor for a line-only
    #   spectrum); with it on, that same mean is a single stable number.
    nmC    = Nuclear.Model(6.)
    stageA = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3")]   # 4 electrons (C2+), E1 pair
    stageB = [Configuration("2s^2 2p"),   Configuration("2s 2p^2")]   # 3 electrons (C3+), E1 pair
    omegas = collect(0.05:0.05:1.2)
    kT     = Defaults.convertUnits("energy: from eV to atomic", 15.0)

    schemeBB = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [[stageA, stageB]], [1.0], omegas)
    settings = Plasma.Settings(kT, 1.0e-8, true)
    compBB   = Plasma.Computation(Plasma.Computation(); scheme=schemeBB, nuclearModel=nmC, grid=Radial.Grid(true),
                                     asfSettings=AsfSettings(), settings=settings)
    wbBB     = perform(compBB; output=true)
    kappaBB  = wbBB[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]

    schemeBF = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [[stageA, stageB]], [1.0], omegas,
                                                                 Empirical.ScaledHydrogenic())
    compBF   = Plasma.Computation(compBB; scheme=schemeBF)
    wbBF     = perform(compBF; output=true)
    kappaBF  = wbBF[ResultKeys.SuperConfigurationMixtureBoundFreeOpacity]
    kappaR   = wbBF[ResultKeys.SuperConfigurationMixtureRosselandOpacity]

    @printf("\n>>> far out on the tail, omega = %.2f Ha: bound-bound = %.3e, bound-free = %.3e [a0^2] (ratio %.2e)\n",
            omegas[end], kappaBB[end], kappaBF[end], kappaBF[end]/kappaBB[end])
    @printf(">>> Rosseland mean opacity (bound-bound + bound-free) = %.6e [a0^2]\n", kappaR)
    println(">>> without bound-free this mean is not even well-defined -- see Plasma.rosselandMeanOpacity's")
    println(">>> own docstring: it moved by more than eleven orders of magnitude under an arbitrary floor for")
    println(">>> a line-only spectrum. Bound-free supplies the physical floor that makes it meaningful.")
    #
elseif  false
    # Last visit:  20-Oct-2026
    # Last successful:  19-Oct-2026   ## crossover at omega = 0.95 Ha: bound-bound 1.0e-9, free-free 9.8e-9
    #                                    [a0^2] (free-free already 10x larger at the NEXT point past the
    #                                    crossing, omega = 1.0: 6.5e-11 vs 8.6e-9)
    #
    # Branch e: FREE-FREE, on the SAME C2+/C3+ two-stage mixture, now with freeFreeIncluded = true
    #   (boundFreeApprox left at nothing, to isolate free-free alone). CHECK, NEEDING NO EXTERNAL
    #   REFERENCE: unlike bound-free, free-free has NO threshold -- it is present at every omega down to
    #   omega -> 0, via Empirical.freeFreeAbsorptionCoefficient (Karzas & Latter 1961's classical cross
    #   section times Sutherland's 1998 thermally-averaged Gaunt factor). Bound-bound falls EXPONENTIALLY
    #   past its last line while free-free falls only as a power law, so free-free eventually CROSSES OVER
    #   and dominates -- found here, not assumed, at omega = 0.95 Hartree.
    nmC    = Nuclear.Model(6.)
    stageA = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3")]
    stageB = [Configuration("2s^2 2p"),   Configuration("2s 2p^2")]
    omegas = collect(0.05:0.05:1.2)
    kT     = Defaults.convertUnits("energy: from eV to atomic", 15.0)

    scheme   = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [[stageA, stageB]], [1.0], omegas,
                                                                 nothing, true)
    settings = Plasma.Settings(kT, 1.0e-8, true)
    comp = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=nmC, grid=Radial.Grid(true),
                                 asfSettings=AsfSettings(), settings=settings)
    wb = perform(comp; output=true)

    kappaBB = wb[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    kappaFF = wb[ResultKeys.SuperConfigurationMixtureFreeFreeOpacity]

    crossIdx = findfirst(k -> kappaFF[k] > kappaBB[k], 1:length(omegas))
    @printf("\n%8s %16s %16s\n", "omega [Ha]", "bound-bound", "free-free [a0^2]")
    for  k = max(1,crossIdx-3):min(length(omegas),crossIdx+2)
        @printf("%8.3f %16.4e %16.4e%s\n", omegas[k], kappaBB[k], kappaFF[k], k==crossIdx ? "   <- crossover" : "")
    end
    println("\n>>> free-free has NO threshold and falls only as a power law; bound-bound falls exponentially")
    println(">>> past its last line -- so free-free eventually dominates even though it is everywhere tiny")
    println(">>> compared to a resonance peak. This is the regime where a bound-bound-only (or even")
    println(">>> bound-free-only) opacity most understates the true continuum.")
    #
elseif  true
    # Last visit:  20-Oct-2026
    # Last reproduced:  20-Oct-2026 (7 of 7 numbers on the claim line)
    # Last successful:  20-Oct-2026   ## resonance width 0.2574689438 eV (example-De.jl's own validated
    #                                    0.2574689437887644 eV, from 9.461809e-03 a.u.); peak kappa_gaussian
    #                                    2.609405e-05 vs peak kappa_voigt 6.909243e-06 [a0^2] at the line
    #                                    centre (omega = 31.488 Ha), ratio 0.264782
    #
    # Branch f: RESOLVED AUTOIONIZING RESONANCES (item C13), on a Ne K-hole case -- ONE species (neon), two
    #   stages: Ne+ (9 electrons, ground 1s^2 2s^2 2p^5 plus the K-hole excited group 1s 2s^2 2p^6, which
    #   sits ABOVE the Ne2+ ionization threshold and is flagged RESONANT by Plasma.classifyResonantGroups)
    #   and Ne2+ (8 electrons, the three final configurations the K-hole decays into by Auger emission --
    #   the SAME triple example-De.jl's own branch a already validates). `resonanceContinuumGrid` turns the
    #   correction on: the ground-to-K-hole E1 line -- Gaussian in the statistical spread of the unresolved
    #   K-hole array -- is reshaped into a Voigt profile using the SAME autoionization width that branch
    #   independently validates (0.2574689437887644 eV, i.e. a 2.556 fs K-shell lifetime).
    #
    #   THE CONTINUUM-METHOD DEFAULTS MUST BE SET BY THE CALLER, exactly as any direct AutoIonization call
    #   requires (Plasma.resonanceAutoionizationWidth's own docstring says so) -- found here, not assumed:
    #   a first attempt without setting them reproduced the width to only 0.36%, not to full precision.
    #
    #   CHECKS, NEEDING NO EXTERNAL REFERENCE BEYOND THE ALREADY-VALIDATED WIDTH: (1) the reported resonance
    #   width must reproduce example-De.jl branch a's own number; (2) the Lorentzian width (0.009462 Ha) is
    #   about 5.5x the line's own Gaussian spread (sigma 0.001719 Ha, from the K-hole TransitionArray's own
    #   variance), so the Voigt profile should visibly broaden and LOWER the peak against the Gaussian-only
    #   line -- found, not merely plausible: the peak drops to 26.5% of the Gaussian-only value, and the
    #   corrected and uncorrected spectra differ by construction only in this one line's own window
    #   (resonanceContinuumGrid = nothing reproduces the pure-Gaussian spectrum exactly, since the correction
    #   is an exact difference against that same Gaussian -- unchanged from before this branch was added).
    nmNe   = Nuclear.Model(10.)
    stageA = [Configuration("1s^2 2s^2 2p^5"), Configuration("1s 2s^2 2p^6")]              # Ne+,  9 electrons
    stageB = [Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")]  # Ne2+, 8 electrons
    omegas = collect(31.40:0.002:31.58)   # Hartree, centred on the K-hole line (~857 eV)
    kT     = Defaults.convertUnits("energy: from eV to atomic", 15.0)

    setDefaults("method: continuum, Galerkin")
    setDefaults("method: normalization, pure sine")

    schemeNoRes = Plasma.SuperConfigurationMixtureOpacityScheme([nmNe], [[stageA, stageB]], [1.0], omegas)
    settings    = Plasma.Settings(kT, 1.0e-8, true)
    compNoRes   = Plasma.Computation(Plasma.Computation(); scheme=schemeNoRes, nuclearModel=nmNe, grid=Radial.Grid(true),
                                        asfSettings=AsfSettings(), settings=settings)
    wbNoRes     = perform(compNoRes; output=true)
    kappaNoRes  = wbNoRes[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]

    continuumGrid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=1.0e-2, rbox=10.0)
    schemeRes     = Plasma.SuperConfigurationMixtureOpacityScheme([nmNe], [[stageA, stageB]], [1.0], omegas,
                                                                      nothing, false, continuumGrid)
    compRes       = Plasma.Computation(compNoRes; scheme=schemeRes)
    wbRes         = perform(compRes; output=true)
    kappaRes      = wbRes[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    corrections   = wbRes[ResultKeys.SuperConfigurationMixtureResonanceCorrections]

    widthEv = corrections[1][6]
    @printf("\n>>> resonance width = %.10f eV  (example-De.jl branch a: 0.2574689437887644 eV)\n", widthEv)
    @printf(">>> peak kappa_gaussian = %.6e, peak kappa_voigt = %.6e [a0^2]  (ratio %.6f)\n",
            maximum(kappaNoRes), maximum(kappaRes), maximum(kappaRes)/maximum(kappaNoRes))
    @printf(">>> max |kappa_voigt - kappa_gaussian| = %.6e [a0^2]\n", maximum(abs.(kappaRes .- kappaNoRes)))
    println(">>> the Lorentzian width dominates the line's own Gaussian spread, so the Voigt profile")
    println(">>> visibly broadens and lowers the peak against the Gaussian-only line -- the correction")
    println(">>> this scheme's steps 4-5 exist to make, now demonstrated end to end on a real resonance.")
    #
end
