# How do I ...?

JAC ships with several hundred worked example branches. This page is **not** an index of them — it is a short,
deliberately incomplete list of the questions users ask most, each pointing at *one* example that answers it.
Start here; the [Examples](examples.md) page is the full working reference when you need something not listed.

Every entry below points at an example branch that carries a `Last successful` date, so the test suite guards it.

## Structure and levels

**How do I get the level structure and energies of an ion?**
`example-Ma.jl` — `computeLevelEnergies(...)` from the `ForPedestrians` module, the shortest path from a
configuration to a level scheme, and the one to start with. `example-Ac.jl` shows the fuller route through a
configuration-interaction expansion for open-shell multiplets, which is what almost everything else on this
page builds on.

**How do I improve the orbitals with correlation?**
`example-Ai.jl` — restricted-active-space (RAS) expansions, adding correlation layer by layer. See also
[Self-consistent fields](scf-routes.md) for which SCF route to use and what it costs.

**How do I choose a radial grid — and know when it is wrong?**
`example-Ab.jl` — B-spline primitives and one-particle spectra in a local potential.
`Basics.recommendedGrid(configs, nuclearModel)` chooses a sensible box for you, and is the recommended way to
build one. A box that is far too *large* starves the basis exactly as badly as one that is too small, and the
symptom looks like an angular-momentum bug rather than a grid problem.

**How do I include the Breit interaction?**
`example-Ad.jl` — the frequency-independent Breit interaction in the CI step.

## Radiative processes

**How do I get transition rates and oscillator strengths (Einstein A)?**
`example-Ca.jl` for the `Einstein` module on a single multiplet, or `example-Da.jl` for `PhotoEmission`
between separately generated initial and final multiplets.

**How do I get a photoexcitation or photoionisation cross section?**
`example-Db.jl` for photoexcitation, `example-Dc.jl` for photoionisation. Both need a continuum orbital, so
the radial grid matters more here than for bound-state work.

**How do I get a radiative-recombination cross section?**
`example-Dd.jl` — `PhotoRecombination`, the inverse of photoionisation.

## Electron-collision processes

**How do I get an Auger rate and a level width?**
`example-De.jl` — `AutoIonization`. The rate is a width; JAC prints both.

**How do I get dielectronic-recombination strengths and rate coefficients?**
`example-Df.jl` — resonant capture followed by radiative stabilisation.

**How do I get an electron-impact excitation cross section?**
`example-Dl.jl` — `ImpactExcitation`.

## Level properties

**How do I get hyperfine A and B constants?**
`example-Cb.jl` — the `Hfs` module, hyperfine parameters and the hyperfine representation.

**How do I get Landé g-factors and Zeeman splittings?**
`example-Cd.jl` — the `LandeZeeman` module.

**How do I get an isotope shift?**
`example-Cc.jl` — mass and field shift from a single multiplet.

## Larger computations

**How do I get charge-state distributions after inner-shell ionisation?**
`example-Fb.jl` — a decay cascade followed through to the final charge states.

**How do I get plasma rate coefficients?**
`example-Jb.jl` — collisional-radiative modelling on top of the atomic data.

## And the question worth asking of every result

**How do I know whether to trust what came back?**

This is the one entry here that is about judgement rather than capability, and JAC answers a good part of it by
itself:

* **Gauge agreement.** Any radiative quantity is computed in both Coulomb and Babushkin gauge. They agree when
  the wave functions are good and diverge when they are not. `PhotoEmission` prints the ratio and a verdict.
* **Cowan's cancellation factor**, printed beside it, says whether the amplitude survived as a small difference
  of large terms — a value below 0.1 means the number carries little information, whatever the gauge agreement.
* **The SCF route reports on itself** — whether it converged, how far it got, and how far its answer can be
  trusted; see [Self-consistent fields](scf-routes.md).

None of these proves a number right: two gauges can agree on a wrong answer when the same correlation is
missing from both. They are reliable in the other direction — a bad indicator is good evidence of a bad number,
and that is what makes them worth reading before you publish anything.
