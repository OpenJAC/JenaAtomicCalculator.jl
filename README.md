

<!-- The codecov badge below deliberately keeps the OLD repository name: codecov keys the project under
     `OpenJAC/JAC.jl` and does not follow the GitHub rename.  Measured 02-Nov-2026, the old path reports 41 %
     while `OpenJAC/JenaAtomicCalculator.jl` reports "unknown".  Do not "correct" it without first checking
     that codecov has been re-pointed. -->
| **Documentation**                                                               | **Build Status**                                                                                | **Try JAC in Cloud** 
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![Documentation](https://img.shields.io/badge/Documentation-stable-blue)](https://openjac.github.io/JenaAtomicCalculator.jl/stable) [![Documentation](https://img.shields.io/badge/Documentation-dev-blue)](https://openjac.github.io/JenaAtomicCalculator.jl/dev) | [![Build Status](https://github.com/OpenJAC/JenaAtomicCalculator.jl/workflows/CI/badge.svg)](https://github.com/OpenJAC/JenaAtomicCalculator.jl/actions) [![codecov](https://codecov.io/gh/OpenJAC/JAC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenJAC/JAC.jl)  | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/OpenJAC/JenaAtomicCalculator.jl/master) |



# Jena Atomic Calculator (JAC) for the computation of atomic representations, processes and cascades


## What is JAC, the JenaAtomicCalculator ?

**JAC**, the **Jena Atomic Calculator**, provides an open-source Julia package for doing atomic computations of various kind
and complexity. In particular, JAC is a (relativistic) electronic structure code for the computation of (atomic many-electron) 
interaction amplitudes, properties as well as a good number of excitation and decay processes for open-shell atoms and ions 
across the whole periodic table. In recent years, moreover, emphasis has been placed on atomic cascade computations
in different physical contexts, on plasma computations (average-atom, ionic mixtures, kinetic and super-configuration
schemes), on symbolic analysis (simplifications) of expressions from Racah's algebra, and on atomic response
computations -- dynamic multipole polarizabilities and dispersion coefficients via a direct inhomogeneous-equation
solve, orchestrated over an atom's active electrons. Further work is ongoing (or still planned) on strong-field
responses to external fields and particles beyond the simplest case, and on the time-evolution of atoms and ions in
the framework of (time-dependent) density matrix and Liouville equations.

A primary guiding philosophy of JAC was to develop a **general and easy-to-use toolbox for the atomic physics community**, 
including an interface that is equally accessible for scientists working in astro and plasma physics, atomic spectroscopy,
theoretical physics as well as for code developers. Beside of its simple use, however, we wish to provide and support a modern
code design, a reasonable detailed documentation of the code and features for integrated testing. Indeed, the JAC toolbox 
facilitates many typical computations and the handling of atomic data by providing input interfaces similar to what one 
uses in a in spoken or written language. Shortly speaking, JAC aims to provide a powerful **platform for daily use and 
to extent atomic theory towards new applications** or, eventually, a **community platform for Just Atomic Computations**.

**Remark**: Although major efforts have been undertaken during the past eight years, JAC is still in an *early state* 
of development and includes features that have been implemented only partly or which are not yet tested well. JAC is an
open-source *physics code* (in its original sense) for which, despite of possible failures and deficiencies, we kindly ask 
potential users and developers for response, support and encouragement.




**Not sure whether JAC already does what you need?** [`CAPABILITIES.md`](CAPABILITIES.md) is a register,
mechanically re-derived from every example branch JAC ships, that lists what JAC can do in plain physics terms,
which module covers it, and the example branch to open first -- including an honest note where something is
not yet verified or not ready for real use. Not a substitute for the documentation below, just a fast first
stop. Search it with `julia tools/searchCapabilities.jl "<term>"`, or just open the file and use your editor's
find.


## *Kinds* of computations

JAC distinguishes and aims to support **different kinds of computations**, summarized below. This list gives
the shape of what JAC does; for the current, itemized detail -- which module, which example branch, and
whether it has been physically verified yet -- see [`CAPABILITIES.md`](CAPABILITIES.md) instead of trying to
keep this list exhaustive.

1. **Atomic computations**, based on explicitly specified electron *configurations*: the frequently applied
    computation of level energies and atomic state representations, together with one or several atomic
    *properties* for selected levels -- hyperfine structure, Landé/Zeeman splitting, isotope shifts, static and
    dynamic multipole polarizabilities and dispersion coefficients, crystal-field and Stark splittings, among
    others -- or **one** selected *process* at a time when levels from two or more multiplets and/or charge
    states are involved, such as photoionization, Auger decay, dielectronic or radiative recombination, and a
    good number of further excitation and decay processes across the periodic table.
2. **Atomic representations**: the generation of different representations of atomic wave functions. This
    includes systematically-enlarged, layer-by-layer restricted active-space (RAS) expansions -- now scaling to
    several thousand CSF per layer, with the self-consistent-field and CI-matrix construction steps
    parallelized -- an explicit second-order (effective-Hamiltonian) treatment that folds a weakly-coupled
    correlation layer directly into the CI matrix instead of enlarging it variationally, approximate Green
    (-function) expansions, and bi-orthogonal transformations between independently-generated orbital sets.
3. **Interactive computations**: the large set of JAC methods can always be used directly at Julia's REPL or
    in short scripts, built from JAC's own data types -- orbitals, configuration-state functions (CSF), atomic
    bases, levels, multiplets, and others -- down to the observable of interest (energies, expansion
    coefficients, amplitudes, rates, cross sections, ...). For a first calculation without assembling these by
    hand, the **`ForPedestrians`** module provides a thin, simplified calling layer (`computeLevelEnergies`,
    `computeTransitionRates`, `computeCrossSections`, `displaySpectrum`, and similar one-call entry points) --
    the fastest on-ramp into JAC, and the natural place to start.
4. **Atomic cascade computations**: an atomic cascade typically connects ions of one element across three or
    more charge states through processes such as photoionization, radiative and dielectronic recombination,
    Auger decay, photoexcitation and emission. A predefined *cascade scheme* selects the state space, the
    processes considered at each step, and the limitations needed to keep the computation feasible; the work is
    split into the cascade computation itself (the many-electron amplitudes) and a separate simulation step that
    combines them for a given experimental scenario -- including, more recently, stepwise-decay, hollow-ion,
    resonant-capture and expansion-opacity schemes.
5. **Atomic response computations**: the response of an atom to an external field or particle beam, beyond a
    single transition amplitude. This now includes a direct inhomogeneous-equation solve for the core
    polarization behind dynamic multipole polarizabilities and dispersion coefficients (`DalgarnoLewis`), an
    orchestration of one such solve per active electron (`RandomPhase`), and a frozen nonlocal (Dirac-Fock)
    exchange operator usable in this context without an SCF iteration. Strong-field ionization in the SFA
    (above-threshold-ionization energy and angular distributions) is supported for hydrogen-like initial
    states; photoelectron momentum distributions in this framework are explicitly postponed -- see
    `CAPABILITIES.md` for which branches are ready to use and which are not.
6. **Plasma computations**: atoms and ions embedded in a plasma environment, where the level structure and
    effective single-electron states remain partly intact. This spans the average-atom model, Saha-Boltzmann
    ionic mixtures in LTE, shifted and broadened lines, dielectronic-satellite diagnostics, kinetic (non-LTE)
    collisional-radiative level populations, super-configuration ionization balance and opacities (including
    multi-species mixtures), and a kinetic charge-balance scheme with a per-process empirical / semi-empirical /
    structure-based rate ladder.
7. **Symbolic evaluation of expressions from Racah's algebra**: the algebraic transformation and
    simplification of Racah expressions -- any number of Wigner n-j symbols together with integrals over
    spherical harmonics, Wigner rotation matrices, and Kronecker/triangular deltas -- based on a large set of
    special values, orthogonality relations and sum rules, cf. the monograph by Varshalovich *et al* (1988).
8. **Empirical computations**: processes such as single and multiple electron-impact ionization, stopping
    powers, tunneling ionization, charge exchange, and several plasma-rate-coefficient models that are not
    efficiently reached by *ab-initio* many-body techniques, computed instead from simple electronic-structure
    input combined with established empirical formulas (module `Empirical`) -- used only where no feasible
    *ab-initio* alternative exists.
9. **Particle- and photon-scattering computations**: elastic electron and positron scattering at atoms (Dirac
    partial waves, differential cross sections, the Sherman function, transport cross sections), one-photon
    positron-electron annihilation, Rayleigh and resonant-inelastic (RIXS) photon scattering, the cheap
    form-factor limit for Rayleigh/Compton scattering, bound-free pair creation, and twisted (Bessel) photon
    beams.
10. **Further, more specialized computations**: muonic atoms (binding energies, muonic X-rays, screening),
    Stark/line-broadening widths via the Dimitrijevic & Konjevic method, three-body recombination rate
    coefficients via detailed balance, deep-learning estimation of missing level energies, and the
    time-evolution of atomic statistical tensors (population and coherences) under the Liouville equation for
    light-pulse-driven dynamics within a fixed level space.


## Installation

JAC is a registered Julia package; install it like any other:
```julia
import Pkg
Pkg.add("JenaAtomicCalculator")
```
and then, in any Julia session:
```julia
using JenaAtomicCalculator
```

JAC supports **Julia 1.10 and 1.12** (see `Project.toml`'s `[compat]` section for the exact bound in force).

Once installed, [`CAPABILITIES.md`](CAPABILITIES.md) and the `ForPedestrians` module (kind 3, above) are the
fastest way to find out whether JAC already does what you need and to get a first calculation running.

## Author and principal developer

JAC is created and maintained by **Stephan Fritzsche** (Helmholtz Institute Jena &
GSI Helmholtz Centre for Heavy Ion Research, Darmstadt).

Contact: s.fritzsche@hi-jena.gsi.de

If you use JAC in published work, please cite the reference article:

> S. Fritzsche, *A fresh computational approach to atomic structures, processes and cascades*,
> Comput. Phys. Commun. **240**, 1–14 (2019).
> DOI: [10.1016/j.cpc.2019.01.012](https://doi.org/10.1016/j.cpc.2019.01.012)

A `CITATION.cff` file is included in this repository; GitHub displays a *Cite this repository*
button in the sidebar that formats the citation automatically.

## How JAC is developed

Most changes to JAC are made with the assistance of an AI coding assistant, working under a written set of
rules kept in [`CLAUDE.md`](CLAUDE.md). [`CONTRIBUTING.md`](CONTRIBUTING.md) states what an external
contribution needs; **["How JAC is developed"](https://openjac.github.io/JenaAtomicCalculator.jl/stable/development/)**
states the AI-assistance record plainly, with the numbers re-derived from `git log` rather than hand-typed.


## Current limitations of JAC

Although JAC has been designed for all atoms and ions across the periodic table, a number of limitations occur:
* Bound-orbital self-consistent-field computations are still based on a local potential (e.g. core-Hartree,
  Kohn-Sham, Dirac-Hartree-Slater, ...) that can be controlled by the user, rather than on a fully nonlocal
  Hartree-Fock-type exchange.
* The CI-matrix construction and the SCF gradient steps were parallelized during 2026, and RAS expansions now
  routinely reach several thousand CSF per layer; very large expansions remain, however, the most expensive
  part of a calculation, and cost still grows quickly with layer count and nuclear charge (see
  `CAPABILITIES.md`'s note on the RAS entry).
* Continuum orbitals are still generated, by default, in a local potential of the ionic core. A frozen
  nonlocal (Dirac-Fock) exchange operator is now available for this case (`Bsplines.setupExchangeMatrix` /
  `Continuum`, used e.g. by the `RandomPhase` response module) without requiring a full SCF iteration -- but it
  is not yet the default everywhere a continuum orbital is generated.

