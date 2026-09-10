# Self-consistent fields: which route, and what to expect

A self-consistent-field computation in JAC answers two separate questions, and it is worth keeping them apart
because they are chosen independently:

* **what the orbitals are optimized FOR** — a mean field, an average level, or one particular level. This is
  the `scField` of [`AsfSettings`](@ref).
* **HOW those orbitals are found** — by which numerical procedure, and with what iteration budget. This is the
  `scfRoute`, and since version 0.6.0 you can choose it.

Most users never need to touch the second one. The default, `Basics.AutomaticRoute()`, takes the standard
procedure for whichever field you asked for, so naming a field remains a complete specification.

## The fields

| `scField` | what it gives you |
|---|---|
| `Basics.DFSField()` | a mean Dirac–Fock–Slater field; the default, and the usual starting point |
| `Basics.ALField()` | an average-level field: one orbital set for a whole multiplet |
| `Basics.EOLField()` | orbitals optimized for one selected level (or a weighted set of levels) |

## The routes

| `scfRoute` | applies to | typical cost |
|---|---|---|
| `Basics.AutomaticRoute()` | any field | the standard choice for that field |
| `Basics.MeanFieldRoute(n)` | DFS, HS, KS, CH | 3–10 iterations |
| `Basics.AverageLevelRoute(n)` | AL | tens of iterations |
| `Basics.RotationRoute(n)` | EOL | hundreds to thousands of iterations |
| `Basics.FockRoute(n)` | EOL | ~10–25 iterations |
| `Basics.NewtonRoute(n)` | — | not implemented; raises and names the alternatives |

Each route carries **its own iteration budget**, because the natural counts differ by two orders of magnitude
between them. Before 0.6.0 a single `AsfSettings.maxIterationsScf` served all of them, which is why its default
of 24 was simultaneously right for a mean field and far too small for an optimized level.

## Choosing an EOL route, and what each is worth

This is the only choice that really requires judgement, and the honest guidance is short.

**Use `RotationRoute` for any number you intend to publish.** It minimizes the level energy directly. It
descends monotonically, keeps the orbital basis orthonormal to about 1e-15, and is the route against which
everything else here was measured. Its cost is iterations: 60–400 on small cases, and more on larger ones.

!!! note "What `CONVERGED` means, and the reference's own uncertainty"
    This route finds a *stationary point*, and which one it finds depends on where it started. Measured on three
    systems from two sensible starting bases, the converged energies differ by **0.05 to 0.6 mHa** — and every run
    reported `CONVERGED`. Asking for a tighter `accuracyScf` does not close the gap: the iteration already stops
    at the energy's own double-precision resolution, so the runs are converged as far as arithmetic allows and
    still differ. The variational principle still applies, so a lower energy is the better calculation, but
    nothing identifies the lowest as global. **Where a milli-Hartree matters, run from more than one starting
    basis and keep the lowest.** The accuracy figures in the table above are differences against this route, so
    they carry this slack too — negligible at Z = 4, about 18 % of the quoted deficit at Z = 26.

**Use `FockRoute` for exploration, scans, and highly charged ions.** It solves the Fock equations to
self-consistency in roughly 10–25 iterations. It is *approximate*: it converges to a slightly different
condition than the energy minimum, and it lands above the rotation route by

<!-- BEGIN generated: route accuracy -->
| system | absolute | relative to the total energy |
|---|---|---|
| Be-like, Z = 4   | 8.9 mHa | 6e-04 |
| Be-like, Z = 10  | 4.4 mHa | 4e-05 |
| Be-like, Z = 26  | 3.2 mHa | 4e-06 |
| Be-like, Z = 92  | 7.1 mHa | 6e-07 |
| C-like, Z = 6    | 10.2 mHa | 3e-04 |

*Measured 10 September 2026 at commit `5976750`; regenerate with `julia --project=. tools/regenerateRouteAccuracy.jl --write`.*
<!-- END generated: route accuracy -->

The pattern is the useful part: **the relative error falls by three orders of magnitude from neutral-like to
highly charged systems.** Near neutrality the competing configurations have comparable weights, the generalized
occupation of a correlating orbital becomes small, and the route can converge onto a poor solution. For a
highly charged ion the reference configuration dominates and the route is excellent.

So: **near-neutral and light — use the rotation route. Highly charged — the Fock route is sound and fast.**

## What the code tells you

You are not expected to remember the table above. The Fock route reports on itself:

* whether it **converged** or exhausted its budget, and at which iteration;
* a standing note that it is fast and approximate, with the measured size of the discrepancy;
* a **specific check that fires by itself** — the two spin-orbit partners of a subshell (say `2p_1/2` and
  `2p_3/2`) must agree in mean radius to within a bound that grows with `Z`, since the splitting is physical.
  A larger deviation is the signature of a failed optimization, and is reported as such.

It also returns the **best** orbital set it saw rather than the last one, since the energy of any orthonormal
set is a variational upper bound.

## Worked example

```julia
using JenaAtomicCalculator

configs = [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]
nModel  = Nuclear.Model(4.0)
grid    = Basics.recommendedGrid(configs, nModel)

# orbitals optimized for the lowest level, by direct minimization -- the accurate route
settings = AsfSettings(AsfSettings();  scField          = Basics.EOLField(),
                                       scfRoute         = Basics.RotationRoute(400),
                                       accuracyScf      = 1.0e-8,
                                       levelSelectionCI = LevelSelection(true, indices=[1]))
multiplet = SelfConsistent.performSCF(configs, nModel, grid, settings)
```

Changing one line selects the fast route instead:

```julia
settings = AsfSettings(settings;  scfRoute = Basics.FockRoute(60))
```

## Restricted active spaces (RAS)

A RAS computation adds correlation layer by layer, and each layer runs an EOL field, so everything above
applies. Two points cost real time if they are overlooked:

**The radial box must be matched to the orbitals.** A box that is far too *large* starves the B-spline basis
exactly as badly as one that is too small — the number of splines is fixed, so a box much wider than the
orbitals spends them on empty space. `Basics.recommendedGrid` chooses a sensible box from the configurations,
and is the recommended way to build one.

**Orbitals from an EOL or RAS basis carry no orbital energy.** A rotation-optimized orbital is not the
eigenfunction of any one-particle operator, so its energy field is `0.0`. Anything that reads an orbital energy
therefore cannot be used on such a basis — frequency-dependent Breit above all, which now refuses explicitly
rather than silently returning the zero-frequency answer.

## Breit interaction

`CoulombBreit(factor)` adds the Breit interaction to the electron–electron interaction. Two practical notes:

* The Breit part is a substantial share of a Coulomb+Breit run — around half, once the Coulomb part has been
  optimized — so budget accordingly.
* `CoulombBreit(factor, :swept)` selects a faster construction of the same quantity, which exploits the fact
  that every Breit kernel factorizes; `:direct` remains the default and the reference form.
