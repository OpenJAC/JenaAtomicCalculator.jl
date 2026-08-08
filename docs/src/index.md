# Jena Atomic Calculator (JAC) for the computation of atomic representations, processes and cascades

*Last update:* August, 8th, 2026

---

**Principal developer:** Stephan Fritzsche —
Helmholtz Institute Jena & GSI Helmholtz Centre for Heavy Ion Research —
[s.fritzsche@hi-jena.gsi.de](mailto:s.fritzsche@hi-jena.gsi.de)

If you use JAC in published work, please cite:
> S. Fritzsche, *A fresh computational approach to atomic structures, processes and cascades*,
> Comput. Phys. Commun. **240**, 1–14 (2019). DOI: [10.1016/j.cpc.2019.01.012](https://doi.org/10.1016/j.cpc.2019.01.012)

---

# What is JAC?

**JAC**, the **Jena Atomic Calculator**, provides an open-source Julia package for doing atomic computations of various kind
and complexity. In particular, JAC is a (relativistic) electronic structure code for the computation of (atomic many-electron) 
interaction amplitudes, properties as well as a good number of excitation and decay processes for open-shell atoms and ions 
across the whole periodic table. In recent years, moreover, emphasis has been placed support atomic cascade computations
in different physical contexts as well as symbolic analysis (simplifications) of expressions from Racah's algebra.
Some further work is done (or planned) to incorporate central features for studying atomic -- strong-field -- responses
to external fields and particles, or the time-evolution of atoms and ions in the framework of (time-dependent) density
matrix and Liouville equations.

A primary guiding philosophy of JAC was to develop a **general and easy-to-use toolbox for the atomic physics community**, 
including an interface that is equally accessible for scientists working in astro and plasma physics, atomic spectroscopy,
theoretical physics as well as for code developers. Beside of its simple use, however, we wish to provide and support a modern
code design, a reasonable detailed documentation of the code and features for integrated testing. Indeed, the JAC toolbox 
facilitates many typical computations and the handling of atomic data by providing input interfaces similar to what one 
uses in a in spoken or written language. Shortly speaking, JAC aims to provide a powerful **platform for daily use and 
to extent atomic theory towards new applications** or, eventually, a **community platform for Just Atomic Computations**.

# *Kinds* of computations

In some more detail, JAC distinguishes and aims to support **different kinds of computations** which can be summarized 
as follows:

1. **Atomic computations**, based on explicitly specified electron *configurations*: This kind refers to the frequently
    applied computation of level energies and atomic state representations. It also supports the computation of either 
    (one or) several atomic properties for selected atomic levels from a given multiplet or of **one** selected process 
    at a time, if atomic levels from two or more multiplets and/or charge states are involved in some atomic transition,
    such as the photo or autoionization of atoms.
2. **Atomic representations**: This kind concerns the generation of different representations of atomic wave functions; 
    in particular, it includes systematically-enlarged restricted active-space (RAS) computations of atomic states 
    and level energies due to a pre-specified active space of orbitals as well as due to the (number and/or kind of) 
    virtual excitations to be taken to be into account. Such RAS computations are normally performed stepwise by making 
    use of the (one-electron) orbital functions from some prior step. Other atomic representations refer to approximate 
    atomic Green functions and, in the future, could be combined with concepts from close-coupling, (exterior) complex 
    scaling, DMRG or perturbation theory.
3. **Interactive computations**: In practice, the (large set of) methods of the JAC toolbox can always be applied also 
    interactively, either directly at Julia's REPL or by using some -- more or less -- short Julia scripts in order to 
    compute and evaluate the desired observables (atomic parameters), such as energies, expansion coefficients, transition
    matrices and amplitudes, rates, cross sections, etc. An interactive computation typically first prepares and applies 
    (certain instances of) JAC’s data types, such as orbitals, configuration-state functions (CSF), atomic bases, levels, 
    multiplets, and others. And like Julia, that is built upon many (high-level) functions and methods, JAC then provides 
    the required language elements for performing specific atomic computations at different degree of complexity and 
    sophistication.
4. **Atomic cascade computations**: An atomic cascade typically includes ions (of one element) in three or more 
    different charge states. These charge states are connected to each other by different atomic processes, such as 
    photo ionization, radiative and dielectronic recombination, Auger decay, the photo excitation and emission, or various 
    others. In practice, moreover, the relative level population of these charge states is usually determined by the
    specific set-up and geometry of the considered experiment. These cascade computations are usually based on some 
    predefined *(cascade) approach* that enables one to automatically select the state-space of the ions, to choose the 
    atomic processes to be considered for the various steps of the cascade, and to specify perhaps additional limitations
    in order to keep the computations feasible. In addition, these cascade computations are generally divided into two
    parts, the (cascade) computation for determining all necessary many-electron amplitudes and the (so-called)
    simulations to combine the amplitudes due to the experimental scenario of interest.
5. **Empirical computations**: Not all atomic  properties and processes, such as the -- single and multiple -- 
    electron-impact ionization, stopping powers or tunnel ionization rates, can be efficiently described by ab-initio 
    many-body techniques. If needed, they are often easier computed by using empirical formulas and models. 
    A number of empirical computations are now supported to deal with such models; are often based on simple
    electronic structure calculations, together with empirically obtained parameters. These computations are only 
    implemented when data are needed but no *ab-initio* computations of the involved processes appears to be feasible. 
6. **Plasma computations**: The notion of atoms and ions in plasma has been frequently applied to analyze the 
    behaviour of plasma in situations where the level structure and effective single-electron states remain partly 
    intact. Useful plasma computations refer to shifted photo lines, ionic mixtures in local and non-local 
    thermodynamic equilibria, or applications of the average-atom model.
7. **Symbolic evaluation of expressions from Racah's algebra**: This kind refers to the algebraic transformation
    and simplification of (Racah) expressions, which may generally include any number of Wigner n-j symbols 
    of different kind as well as (various integrals over) the spherical harmonics, the Wigner rotation matrices
    and the Kronecker and triangular deltas. Of course, the complexity of such *Racah expressions* increases 
    very rapidly as more Wigner symbols are involved. A symbolic evaluation of these expressions is naturally 
    based on the knowledge of a large set of special values, orthogonality relations and *sum rules* that may include 
    rules with a (multiple) summation over dummy indices, cf. the monography by Varshalovich *et al* (1988).

    
Several other *different kinds of computations* have been prepared and will support the applications of the 
JAC toolbox but come with a rather limited implementation so far.
    
8. **Atomic responses**: With this kind, we partly support computations in intense laser field; they also help 
    analyze the response of atoms to incident beams of light pulses and particles, such as field-induced 
    ionization processes, high-harmonic generation and several others. For these responses, the detailed 
    structure of atoms and ions has not been considered much until today. A partial-wave formulation of these 
    strong-field processes enables one to clearly distinguish between contributions due to the atomic target,
    the Volkov states, or the shape and phase of the incident light.
9. **Atomic time-evolution of statistical tensors**: We here wish to simulate the population and coherences
    of (atomic) levels using the *Liouville equation*, when atoms and ions are irradiated by (intense) light
    pulses. For these computations, however, we shall assume always that the level structure of the atoms is kept 
    intact. Further (decay) processes of the excited atoms and ions can be taken into account by some *loss 
    rate*, but without that the atoms can leave the *pre-specified space of sublevels*. In particular, we here 
    plan to consider the interaction of atoms and ions with pulses of different shape, polarization strength 
    and duration.


# Installation

The development version of JAC can installed installed directly in the interactive Julia REPL (Read-Eval-Print-Loop)
Press `]` to enter the Julia package mode and type `add https://github.com/OpenJAC/JAC.jl` to install JAC.

```
julia> ]
(@v1.10) pkg> add https://github.com/OpenJAC/JAC.jl
```

# Development

As a **Scientific package** users are encouraged to explore the capabilities of JAC, and modify the package as per their 
requirement. JAC comes with an MIT 'Expat' [License](@ref jac-license).

For development purpose, user can install JAC from GitHub to their **desired direcory** using `git clone` in the 
terminal. This will download a copy of JAC from the development branch of GitHub repository.
```
git clone https://github.com/OpenJAC/JAC.jl.git
```

Next `cd` to the `JAC.jl` directory and start a new Julia session

```
shell> cd JAC.jl
shell> julia
```

Then activate the `JAC.jl` environment and install the dependancy packages of `JAC.jl`
```
julia> ]                   # Change to the Julia package mode
(@v1.10) pkg> activate .
(JAC) pkg> instantiate     # Installs the dependancy packages of JAC
(JAC) pkg> resolve         # Required **Only if** there is a cinflict in the dependancies
(JAC) pkg> develop .       # Adds the current directory path to Julia local Registry
```

# User Guide and Manual

A detailed 
[User Guide, Compendium & Theoretical Background to JAC](https://github.com/AlokaSahoo/JAC.jl/blob/master/docs/UserGuide-Jac.pdf)  
is available that describes the **use and underlying atomic theory** of the JAC code.


# Dependencies 

Dependencies and external code used in JAC

The JAC code makes use of:
* standard Julia packages, such as BSplineKit, SpecialFunctions, GaussQuadrature, GSL and QuadGK.
* Matrix elements from G. Gaigalas and S. Fritzsche, Comp. Phys. Commun. 267, 108086 (2021).


# Current limitations of JAC

Although JAC has been designed for all atoms and ions across the periodic table, a number of limitations occur:
* All self-consistent-field computations are based so far on a local potential (e.g. core-Hartree, Kohn-Sham, 
  Dirac-Hartree-Slater, ...) that can be controlled by the user.
* Until the present, no serious optimization has been done for the code; this restricts most computations
  to CSF expansion with several hundred CSF.
* All continuum orbitals are generated in a local potential (Dirac-Hartree-Slater or others, as above) of the ionic core, 
  and without the explicit treatment of the exchange interaction.

