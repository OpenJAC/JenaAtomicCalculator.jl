# JAC capabilities register

**What this file is for.** JAC has grown past the point where any one person -- maintainer, collaborator, or an AI tool reading the repository -- can hold its full feature set in their head. This file exists to answer two questions quickly, for a human or for an AI assistant acting on a human's behalf:

1. **Is JAC already able to do what I need, even approximately?**
2. **Which one or two example branches should I open first to try it and adapt it?**

**What this file is NOT.** It is not the manual, and completeness of DETAIL is explicitly not the goal -- simple use and fast orientation are. The full technical detail lives in each module's own docstrings and in `docs/`; this file only tells you whether to look there at all, and where to start. It IS complete in COVERAGE, though: every branch listed by `examples/examples.jl` gets its own row here, mechanically, so nothing in that index goes missing from this one.

**How to trust what it says.** Every row's date comes from that branch's own `# Last successful: DD-Mon-YYYY` comment line in `examples/*.jl` -- the one marker in this codebase that a human, not a tool, has ever written after physically judging the numbers (see `CLAUDE.md`, Rule 7). A row marked **unknown/undated** means the module exists and the branch runs, but nobody has yet judged the physics, or an earlier judgement was withdrawn -- **undated is not the same as unimplemented or abandoned, and the text after the dash says which it is wherever the branch's own author already wrote one** (item 27, 17-Nov-2026): a branch deliberately parked, one that ran and failed, a guard test that is SUPPOSED to not pass, and one nobody has simply gotten to yet all now read differently, taken verbatim from that branch's own `# Last successful: unknown ...` comment rather than flattened to one generic sentence. Where no such note exists in the source, the generic sentence is all there honestly is to say. An omitted `Module:` line means the same kind of honesty -- this register's generator could not confidently name one, rather than guessing.

**Restrictions column.** A one-line hint, not a specification. It exists to stop a false start, not to replace reading the branch. A blank entry means no restriction is currently flagged -- not that none exists.

**One file, several branches.** Several example files carry more than one `if true`/`if false` branch, each dated independently. The date shown below is the MOST RECENT judged date found ANYWHERE in that file, which may belong to a different branch than the one you actually want -- always check the specific branch's own `Last successful` comment before trusting a copied date. `example-Ar.jl` is a real case of this: one branch there carries a genuine 13-Sep-2026 date while a different, still-undated branch in the SAME file is the one most readers would reach for first.

**Searching this file.** Use your editor's find, a plain `grep -i -B1 -A5 "<term>" CAPABILITIES.md`, or the small helper that prints a matched entry with its full context: `julia tools/searchCapabilities.jl "<term>"`. All three are plain substring search -- "DR" will also match "quaDRupole"; that is a known, accepted limit of keeping the tool genuinely simple, not a defect.

**Once you have found your branch, the next question is how to turn it into a number worth quoting** -- and that is a different document. `docs/UserSteps-JAC.txt` walks an occasional user through ten steps from a physics question to a publishable result: classifying the problem against JAC's seven kinds of computation, choosing the mean field and the correlation model as a PHYSICS decision rather than a default, running the smallest case first and then scaling ONE axis at a time, reading JAC's own printed verdicts before reading the number itself, and finally comparing against something JAC did not compute. This register answers WHETHER and WHERE; those ten steps answer HOW -- and their step 2 is the one that sends a reader back here.

---

## A) Examples, tests & development of the electronic structure part.

### Ab) Apply & test the  a SCF field: B-spline primitives and one-particle spectra in a local potential.

- Start here: `examples/example-Ab.jl`, dated **29-Jul-2026**

### Ac) Apply & test the CI part for open-shell multiplets of increasing complexity, without Breit interaction.

- Start here: `examples/example-Ac.jl`, dated **31-Aug-2026**

### Ad) Apply & test for the frequency-independent Breit interaction for an internally generated neon multiplet.

- Start here: `examples/example-Ad.jl`, dated **15-Aug-2026**

### Ae) Apply & test the QED model corrections to the level structure of atoms and ions.

- Start here: `examples/example-Ae.jl`, dated **29-Jul-2026**

### Af) Generate, normalize & test for continuum orbitals in a local potential.

- Start here: `examples/example-Af.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Ag) Apply & test the jj-LS transformation of levels from a given multiplet.

- Start here: `examples/example-Ag.jl`, dated **18-Aug-2026**

### Ah) Apply & test a mean-field basis, mean-field multiplets as well as a configuration-interaction (CI) expansions.

- Start here: `examples/example-Ah.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Ai) Apply & test for restricted-active-space (RAS) expansions.

- Module: `AtomicState.RasExpansion`
- Start here: `examples/example-Ai.jl`, dated **17-Oct-2026**
- Keywords: RAS, restricted active space, correlation layer
- Restriction: cost grows fast with layer count and nuclear charge -- see `docs/` or ask before a 4th layer on a heavy ion.

### Aj) Apply & test for a Green (-function) expansion.

- Start here: `examples/example-Aj.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Ak) Apply & test the computation of spin-angular coefficients.

- Start here: `examples/example-Ak.jl`, dated **30-Aug-2026**

### Am) Apply & test procedures to establish a configuration-based language and to deal with electron configurations.

- Start here: `examples/example-Am.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### An) Apply & test the bi-orthogonal transformation of two independently-generated multiplets.

- Start here: `examples/example-An.jl`, dated **09-Nov-2026**

### Ao) Apply & test the average-level (AL) and optimized-level (EOL) self-consistent fields, and where they differ.

- Module: `Basics (AverageLevelRoute / RotationRoute -- the AL/EOL SCF routes)`
- Start here: `examples/example-Ao.jl`, dated **01-Sep-2026**

### Ap) Apply & test the three-parameter Fermi nucleus: the shape of the nuclear charge distribution.

- Start here: `examples/example-Ap.jl`, dated **12-Aug-2026**

### Aq) Compare the spin-angular coefficients of SpinAngular against SpinAngular and against GRASP2018.

- Start here: `examples/example-Aq.jl`, dated **25-Aug-2026**

### Ar) Apply & test the second-order treatment of the Q space of a RAS step.

- Module: `Hamiltonian.performCIwithSecondOrderQ`
- Start here: `examples/example-Ar.jl`, dated **13-Sep-2026**

### As) Apply & test the DalgarnoLewis module: a direct inhomogeneous-equation solve in place of an explicit sum over intermediate states.

- Module: `DalgarnoLewis`
- Start here: `examples/example-As.jl`, dated **14-Oct-2026**

### At) Apply & test the RandomPhase module: orchestrating one DalgarnoLewis solve per active electron.

- Module: `RandomPhase`
- Start here: `examples/example-At.jl`, dated **16-Oct-2026**

### Au) Apply & test that the SCF+RAS multi-threading (steps 1-6) changes no result on a mid-size EOL case, only wall time.

- Start here: `examples/example-Au.jl`, dated **26-Oct-2026**

---

## B) Examples, tests & development of atomic amplitudes

### Ba) Apply & test the dipole, em multipole and momentum-transfer amplitudes.

- Module: `Basics (dipole / multipole amplitude routines)`
- Start here: `examples/example-Ba.jl`, dated **28-Aug-2026**

### Bb) Apply & test the parity non-conservation, Schiff moment and anapole moment amplitudes.

- Start here: `examples/example-Bb.jl`, dated **28-Aug-2026**

---

## C) Examples, tests & development of atomic properties

### Ca) Apply & test the Einstein module with ASF from an internally generated multiplet.

- Module: `Einstein`
- Start here: `examples/example-Ca.jl`, dated **30-Aug-2026**

### Cb) Apply & test the Hfs module for HFS A,B parameters and hyperfine representation with ASF from an internally generated multiplet.

- Module: `Hfs`
- Start here: `examples/example-Cb.jl`, dated **28-Aug-2026**
- Keywords: HFS, hyperfine, A constant, B constant
- Restriction: a subshell-factor convention fix landed in 2026 (`CLAUDE.md` Rule 18) -- a number computed before that fix is not comparable to one after it.

### Cc) Apply & test the IsotopeShift module with ASF from an internally generated multiplet.

- Module: `IsotopeShift`
- Start here: `examples/example-Cc.jl`, dated **13-Oct-2026**

### Cd) Apply & test the LandeZeeman module with ASF from an internally generated multiplet.

- Module: `LandeZeeman`
- Start here: `examples/example-Cd.jl`, dated **28-Aug-2026**

### Ce) Apply & test the FormFactor module with ASF from an internally generated multiplet.

- Module: `FormFactor`
- Start here: `examples/example-Ce.jl`, dated **31-Jul-2026**

### Cf) Apply & test the reduced 1-particle density matrix, natural orbitals and related quantities.

- Start here: `examples/example-Cf.jl`, dated **23-Jul-2026**

### Cg) Apply & test the DecayYield module with ASF from an internally generated multiplet.

- Module: `DecayYield`
- Start here: `examples/example-Cg.jl`, dated **13-Oct-2026**

### Ci) Apply & test the AlphaVariation module with ASF from an internally generated multiplet.

- Module: `AlphaVariation`
- Start here: `examples/example-Ci.jl`, dated **21-Oct-2026**

### Cj) Apply & test the MultipolePolarizibility module with ASF from an internally generated multiplet.

- Module: `MultipolePolarizibility`
- Start here: `examples/example-Cj.jl`, dated **14-Oct-2026**

### Ck) Apply & test the CrystalField module (point-charge Stark splitting) with ASF from an internally generated multiplet.

- Module: `CrystalField`
- Start here: `examples/example-Ck.jl`, dated **28-Aug-2026**

### Cl) Apply & test the StarkShift module with ASF from an internally generated multiplet.

- Module: `StarkShift`
- Start here: `examples/example-Cl.jl`
- Date: **unknown/undated** -- RUNS (exit 0) and is the CLEAREST demonstration of why this file

### Cm) Apply & test the StarkZeeman module (exact ASF-basis diagonalization of level mixing in static E and/or B fields of arbitrary relative direction).

- Module: `StarkZeeman`
- Start here: `examples/example-Cm.jl`, dated **01-Aug-2026**

### Cn) Apply & test the WeakInteractionEnhancement module: the parity-non-conserving E1 amplitude and the electron-EDM enhancement factor, both as sums over caller-supplied intermediate states.

- Module: `WeakInteractionEnhancement`
- Start here: `examples/example-Cn.jl`, dated **05-Sep-2026**

### Co) Statistical tensors of an atomic ensemble: alignment and orientation from photoexcitation, the axis they are stated about, and the tensor machinery itself.

- Start here: `examples/example-Co.jl`, dated **28-Aug-2026**

### Cp) Apply & test the NuclearShielding module: the quadrupole antishielding factor of a closed-shell ion, i.e. how its own electrons amplify at the nucleus a field gradient applied from outside.

- Module: `NuclearShielding`
- Start here: `examples/example-Cp.jl`, dated **06-Oct-2026**

### Cq) Multipole polarizabilities alpha^(Mp,K) at real and imaginary frequency, and the long-range dispersion coefficients C_6, C_8, C_10 and the atom-wall C_3 built from them.

- Start here: `examples/example-Cq.jl`, dated **23-Oct-2026**

### Cr) Apply & test the HfsPaschenBack module (exact Hfs.HfLevel-basis diagonalization of hyperfine Paschen-Back/Breit-Rabi mixing in a static B and/or E field).

- Module: `HfsPaschenBack`
- Start here: `examples/example-Cr.jl`, dated **26-Oct-2026**
- Keywords: HFS, Breit-Rabi, Paschen-Back

---

## D) Examples, tests & development of basic atomic processes

### Da) Apply & test the PhotoEmission module with ASF from an internally generated initial- and final-state multiplet.

- Module: `PhotoEmission`
- Start here: `examples/example-Da.jl`, dated **28-Aug-2026**

### Db) Apply & test the PhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.

- Module: `PhotoExcitation`
- Start here: `examples/example-Db.jl`, dated **28-Aug-2026**

### Dc) Apply & test the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.

- Module: `PhotoIonization`
- Start here: `examples/example-Dc.jl`, dated **02-Aug-2026**

### Dd) Apply & test the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.

- Module: `PhotoRecombination`
- Start here: `examples/example-Dd.jl`, dated **28-Aug-2026**

### De) Apply & test the AutoIonization module with ASF from an internally generated initial- and final-state multiplet.

- Module: `AutoIonization`
- Start here: `examples/example-De.jl`, dated **17-Oct-2026**
- Keywords: Auger
- Restriction: a modelling edge case for this process is under active review; check for a more recent date before relying on a result close to a known difficult case (e.g. near-degenerate continua).

### Df) Apply & test the DielectronicRecombination module with ASF from internally generated initial-, intermediate and final-state multiplets.

- Module: `DielectronicRecombination`
- Start here: `examples/example-Df.jl`, dated **01-Sep-2026**
- Keywords: DR
- Restriction: this file carries several branches that are still undated -- check the date of the SPECIFIC branch you copy, not just the file.

### Dg) Apply & test the RayleighCompton module with ASF from an internally generated initial- and final-state multiplets.

- Module: `RayleighCompton`
- Start here: `examples/example-Dg.jl`
- Date: **unknown/undated** -- BLOCKED by a defect in src/module-RayleighCompton.jl itself, see below.

### Dh) Apply & test the MultiPhotonTransition module: TWO-photon emission, absorption and bichromatic absorption.

- Module: `MultiPhotonTransition`
- Start here: `examples/example-Dh.jl`, dated **01-Nov-2026**

### Du) Apply & test the MultiPhotonTransition module: THREE-photon processes (emission sharings; elementary absorption).

- Module: `MultiPhotonTransition`
- Start here: `examples/example-Du.jl`, dated **08-Aug-2026**

### Di) Apply & test the DoubleAutoIonization module with ASF from an internally generated initial- and final-state multiplet.

- Module: `DoubleAutoIonization`
- Start here: `examples/example-Di.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Dj) Apply & test the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.

- Module: `PhotoDoubleIonization`
- Start here: `examples/example-Dj.jl`
- Date: **unknown/undated** -- and it stays unknown: a Coulomb/Babushkin ratio of 1.76 is not a

### Dk) Apply & test the RadiativeAuger module with ASF from an internally generated initial- and final-state multiplet.

- Module: `RadiativeAuger`
- Start here: `examples/example-Dk.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Dl) Apply & test the ImpactExcitation module with ASF from an internally generated initial- and final-state multiplet.

- Module: `ImpactExcitation`
- Start here: `examples/example-Dl.jl`, dated **13-Aug-2026**
- Restriction: at least one tested case has come in well below the high-energy Bethe benchmark; treat an absolute rate as a first estimate until checked against an independent source.

### Dm) Apply & test the InternalRecombination module with ASF from an internally generated initial and final-state multiplet.

- Module: `InternalRecombination`
- Start here: `examples/example-Dm.jl`
- Date: **unknown/undated** -- (see REPORT -- root cause #2 above still open, absolute magnitude not

### Dn) Apply & test the TwoElectronOnePhoton (TEOP) module: both computational strategies compared for the same transition.

- Start here: `examples/example-Dn.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Do) Apply & test the HyperfineInduced module: electronic quenching, nuclear hyperfine mixing, and mixed cases.

- Module: `HyperfineInduced`
- Start here: `examples/example-Do.jl`, dated **06-Aug-2026**
- Keywords: HFS

### Dp) Coulomb IONIZATION of ions by fast ion impact: energy-differential cross sections and the alignment of the residual ion.

- Start here: `examples/example-Dp.jl`, dated **25-Aug-2026**

### Dq) Apply & test the InternalConversion module with ASF from an internally generated initial- and final-state multiplet.

- Module: `InternalConversion`
- Start here: `examples/example-Dq.jl`, dated **15-Aug-2026**

### Dr) Apply & test the CoulombIonization module with ASF from an internally generated initial- and final-state multiplet.

- Module: `CoulombIonization`
- Start here: `examples/example-Dr.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Ds) Apply & test the CoulombExcitation module with ASF from an internally generated initial- and final-state multiplet.

- Module: `CoulombExcitation`
- Start here: `examples/example-Ds.jl`, dated **28-Aug-2026**

### Dt) Apply & test the CrystalFieldEmission module (crystal-field-resolved transitions) with initial- and final-state multiplets declared directly via Atomic.Computation.

- Module: `CrystalFieldEmission`
- Start here: `examples/example-Dt.jl`, dated **28-Aug-2026**

### Dw) Apply & test the GeneralizedOscillatorStrength module: generalized oscillator strengths f_n(K) of the Bethe theory, for bound-bound transitions.

- Module: `GeneralizedOscillatorStrength`
- Start here: `examples/example-Dw.jl`, dated **28-Aug-2026**

### Dx) Apply & test the PhotoRecombinationInterference module: RR and DR amplitudes added coherently; interference in the cross section, beta_2 and the photon polarization.

- Module: `PhotoRecombinationInterference`
- Start here: `examples/example-Dx.jl`, dated **20-Aug-2026**
- Keywords: DR, RR

### Dy) Apply & test the PhotoIonizationInterference module: direct PI and resonant (excitation+autoionization) amplitudes added coherently; the Fano interference term in the cross section.

- Module: `PhotoIonizationInterference`
- Start here: `examples/example-Dy.jl`, dated **20-Oct-2026**
- Keywords: Fano

### Dz) Apply & test the PhotoIonizationPseudostates module: a diagonalize-once CI mixing bound autoionizing CSFs with a box-discretized continuum pseudostate ladder -- the alternative to Dy's amplitude-sum route.

- Module: `PhotoIonizationPseudostates`
- Start here: `examples/example-Dz.jl`, dated **05-Nov-2026**

---

## E) Examples, tests & development of composed atomic processes

### Ea) Apply & test the BeamPhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.

- Module: `BeamPhotoExcitation`
- Start here: `examples/example-Ea.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Eb) Apply & test the  PhotoExcitationFluores module with ASF from an internally generated initial-, intermediate and final-state multiplets.

- Module: `PhotoExcitationFluores`
- Start here: `examples/example-Eb.jl`
- Date: **unknown/undated** -- PARKED under Rule 13 (31-Aug-2026): PARKED BY THE MAINTAINER.

### Ec) Apply & test the PhotoExcitationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.

- Module: `PhotoExcitationAutoIon`
- Start here: `examples/example-Ec.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Ed) Apply & test the Photoionization and PlasmaShift module with ASF from an internally generated initial- and final-state multiplet.

- Start here: `examples/example-Ed.jl`, dated **09-Nov-2026**

### Ef) RETIRED -- plasma-shifted Auger rates moved to Plasma.LineShiftScheme; see example-Jb.jl branch 1.

- **RETIRED** -- kept here only so a search does not come up empty; follow the pointer above.

### Ei) Two-photon IONIZATION of a one-electron system: generalized cross sections in GM, the convergence of the intermediate sum, and the resonance guard.

- Start here: `examples/example-Ei.jl`, dated **01-Sep-2026**

### Ek) Test of the PairAnnihilation1Photon module with ASF from an internally generated initial- and final-state multiplet.

- Module: `PairAnnihilation1Photon`
- Start here: `examples/example-Ek.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Em) Test of the PhotoIonizationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.

- Module: `PhotoIonizationAutoIon`
- No dedicated `examples/example-Em.jl` file was found for this branch code.

---

## F) Examples, tests & development of atomic cascades

### Fa) Stepwise decay cascades: a fast Mg K-hole reference case, and the larger Si^- 1s-3p case.

- Start here: `examples/example-Fa.jl`, dated **23-Aug-2026**
- Keywords: cascade, Auger cascade

### Fb) Cascade SIMULATIONS: one computation, several properties derived from the same cascade data.

- Start here: `examples/example-Fb.jl`, dated **23-Aug-2026**

### Fc) Cascade.PhotoExcitationScheme: photo-excitation cascades and the configuration pre-filter.

- Module: `Cascade.PhotoExcitationScheme`
- Start here: `examples/example-Fc.jl`, dated **06-Aug-2026**

### Fd) Cascade.PhotoIonizationScheme: photo-ionization cross sections of Ne^+, calibrated on neutral Ne.

- Module: `Cascade.PhotoIonizationScheme`
- Start here: `examples/example-Fd.jl`, dated **24-Aug-2026**
- Keywords: cascade

### Fe) Cascade.PhotoAbsorptionScheme: photoabsorption of C^+ near the 1s -> 2p resonances.

- Module: `Cascade.PhotoAbsorptionScheme`
- Start here: `examples/example-Fe.jl`, dated **23-Aug-2026**

### Ff) Cascade.DielectronicRecombinationScheme: KLL dielectronic recombination of helium-like carbon.

- Module: `Cascade.DielectronicRecombinationScheme`
- Start here: `examples/example-Ff.jl`, dated **06-Aug-2026**
- Keywords: cascade, KLL

### Fg) Cascade.RadiativeRecombinationScheme: radiative recombination of helium-like carbon.

- Module: `Cascade.RadiativeRecombinationScheme`
- Start here: `examples/example-Fg.jl`, dated **23-Aug-2026**

### Fh) Cascade.ImpactExcitationScheme: electron-impact excitation of hydrogen-like carbon, direct and resonant.

- Module: `Cascade.ImpactExcitationScheme`
- Start here: `examples/example-Fh.jl`, dated **18-Aug-2026**

### Fi) Cascade.ElectronIonizationScheme: excitation-autoionization of lithium-like carbon.

- Module: `Cascade.ElectronIonizationScheme`
- Start here: `examples/example-Fi.jl`, dated **24-Aug-2026**

### Fj) Cascade.HollowIonScheme: formation and decay of a hollow carbon ion.

- Module: `Cascade.HollowIonScheme`
- Start here: `examples/example-Fj.jl`, dated **23-Aug-2026**

### Fk) Cascade.ExpansionOpacityScheme: bound-bound expansion opacities of Sr^+ in kilonova ejecta.

- Module: `Cascade.ExpansionOpacityScheme`
- Start here: `examples/example-Fk.jl`, dated **18-Aug-2026**

### Fl) Cascade.ElectronIonizationScheme, resonant channels: resonant electron capture with sequential or simultaneous double autoionization, for lithium-like carbon.

- Module: `Cascade.ElectronIonizationScheme`
- Start here: `examples/example-Fl.jl`, dated **24-Aug-2026**

### Fm) Cascade.DielectronicCaptureScheme: resonant capture of a free electron into doubly-excited levels, for the KLL group of helium-like carbon.

- Module: `Cascade.DielectronicCaptureScheme`
- Start here: `examples/example-Fm.jl`, dated **23-Aug-2026**
- Keywords: cascade, KLL

---

## G) Examples, tests & development of symbolic evaluations of Racah expressions

### Ga) Symbolic evaluation by means of special values & recursion relations of the Wigner n-j symbols.

- Start here: `examples/example-Ga.jl`, dated **13-May-2024**

### Gb) Symbolic simplification of Racah expressions by means of sum rules.

- Start here: `examples/example-Gb.jl`, dated **13-May-2024**

### Gc) Symbolic simplification of recoupling coefficients by means of sum rules.

- Start here: `examples/example-Gc.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Gd) Symbolic simplification of Racah expressions including spherical harmonics and Wigner rotation matrices.

- Start here: `examples/example-Gd.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Ge) Symbolic simplification of spherical tensor operators, matrix elements and spherical amplitudes.

- Start here: `examples/example-Ge.jl`, dated **01-Sep-2026**

---

## H) Examples, tests & development of semiempirical computations

### Hb) Compute & test (empirical) impact-ionization cross sections.

- Start here: `examples/example-Hb.jl`, dated **10-Aug-2026**

### Hc) Apply & test semiempirical binding energies, Slater Zeff, Rydberg series, and NIST ionization potentials (module Semiempirical).

- Module: `Semiempirical`
- Start here: `examples/example-Hc.jl`, dated **15-Jul-2026**
- Keywords: Zeff, Rydberg, quick estimate

### Hd) Apply & test the K-shell fluorescence and Auger yields (KrauseAdopted2016) and the empirical Auger rates.

- Start here: `examples/example-Hd.jl`, dated **16-Jul-2026**

---

## I) Examples, tests & development of atomic response computations

### Ib) Apply & test the StrongField module to calculate ATI energy distributions in the SFA with hydrogen-like initial states.

- Module: `StrongField`
- Start here: `examples/example-Ib.jl`
- Date: **unknown/undated** -- PARKED under Rule 13 (26-Aug-2026): POSTPONED BY THE MAINTAINER UNTIL 2027, together with examples/example-Id.
- Keywords: ATI, strong field, SFA
- Restriction: shares the StrongField module with Id; check this branch's own date independently before trusting it.

### Ic) Apply & test the StrongField module to demonstrate the Coulomb asymmetry in ATI azimuthal angular distributions.

- Module: `StrongField`
- Start here: `examples/example-Ic.jl`
- Date: **unknown/undated** -- PARKED under Rule 13 (26-Aug-2026): POSTPONED BY THE MAINTAINER UNTIL 2027, together with examples/example-Id.
- Keywords: ATI, strong field, SFA
- Restriction: shares the StrongField module with Id; check this branch's own date independently before trusting it.

### Id) Apply & test the StrongField module to compute in ATI photoelectron momentum distributions.   [POSTPONED until 2027]

- Module: `StrongField`
- Start here: `examples/example-Id.jl`
- Date: **unknown/undated** -- PARKED under Rule 13 (26-Aug-2026): POSTPONED BY THE MAINTAINER UNTIL 2027, together with examples/example-Id.
- Keywords: ATI, strong field, SFA
- Restriction: explicitly marked [POSTPONED until 2027] in examples.jl itself -- do not use for a real calculation.

---

## J) Examples, tests & development of plasma computations

### Ja) Apply & test the average-atom computations.

- Start here: `examples/example-Ja.jl`, dated **12-Oct-2026**
- Keywords: average atom, plasma, finite temperature
- Restriction: the free-free contribution near zero frequency is an open point; do not rely on this model there without checking for a newer date.

### Jb) Apply & test the line-shift computations.

- Start here: `examples/example-Jb.jl`, dated **28-Aug-2026**

### Jc) Apply & test the Saha-Boltzmann computations for an ionic mixture in LTE.

- Start here: `examples/example-Jc.jl`, dated **23-Oct-2026**
- Keywords: Saha-Boltzmann, LTE, ionization balance

### Jd) Apply & test the Plasma.SatelliteDiagnosticScheme (dielectronic-satellite-to-parent-line intensity-ratio Te diagnostic).

- Module: `Plasma.SatelliteDiagnosticScheme`
- Start here: `examples/example-Jd.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Je) Apply & test the Plasma.CollisionalRadiativeScheme: kinetic (non-LTE) level populations for a small level set.

- Module: `Plasma.CollisionalRadiativeScheme`
- Start here: `examples/example-Je.jl`, dated **04-Aug-2026**

### Jf) Apply & test the Plasma.SuperConfigurationOpacityScheme (populations, oscillator strengths, and a summed opacity spectrum).

- Module: `Plasma.SuperConfigurationOpacityScheme`
- Start here: `examples/example-Jf.jl`, dated **20-Oct-2026**

### Jg) Apply & test the Plasma.SuperConfigurationIonizationScheme (equilibrium ionization balance across configurations).

- Module: `Plasma.SuperConfigurationIonizationScheme`
- Start here: `examples/example-Jg.jl`, dated **17-Oct-2026**

### Jh) Apply & test the Plasma.SuperConfigurationMixtureOpacityScheme (a multi-species mixture, Planck mean opacity).

- Module: `Plasma.SuperConfigurationMixtureOpacityScheme`
- Start here: `examples/example-Jh.jl`, dated **20-Oct-2026**

### Ji) Apply & test the Plasma.KineticChargeBalanceScheme (kinetic population balance across several charge states, with a per-process empirical/semi-empirical/structure-based rate-treatment ladder).

- Module: `Plasma.KineticChargeBalanceScheme`
- Start here: `examples/example-Ji.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

---

## K) Examples, tests & development of deep-learning applications

### Ka) Train & apply a neural network for estimating missing energies in the NIST tables.

- Start here: `examples/example-Ka.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Kb) Apply & test the AtomicFeatures module: the FEATURE VECTORS a neural network is trained on.

- Module: `AtomicFeatures`
- Start here: `examples/example-Kb.jl`, dated **01-Sep-2026**

---

## L) Examples, tests & development of Liouville (time-evolution) computations

### La) Apply & test the Liouville (time-evolution) computations for the stimulated Raman scattering.

- Start here: `examples/example-La.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

---

## M) Examples, tests & development of the ForPedestrians module

### Ma) Test  computeLevelEnergies(ForGivenConfigs/ForIsoelectronicSequence(), ...)  for several atoms and sequences.

- Start here: `examples/example-Ma.jl`, dated **27-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Mb) Test  computeTransitionRates(),  computeBranchingFractions()  and  displaySpectrum()  for ForPhotoEmission and ForAutoIonization.

- Start here: `examples/example-Mb.jl`, dated **28-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Mc) Test  computeCrossSections(ForPhotoIonization(), ...)  and  computeCrossSections(ForPhotoRecombination(), ...)  for several atoms.

- Start here: `examples/example-Mc.jl`, dated **26-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Md) Test  computeLifetimes(ForPhotoEmission(), ...)  and  computeLifetimes(ForAutoIonization(), ...)  for Ne and Ar K-holes; fluorescence yield.

- Start here: `examples/example-Md.jl`, dated **26-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Me) Test  displayCouplings(FineStructure(), ...)  and  displayCouplings(FineStructureLS(), ...)  for Ar^2+ with two inner-shell holes.

- Start here: `examples/example-Me.jl`, dated **26-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Mf) Test  estimateCrossSections(ForImpactIonization(), ...)  for neon and argon; relativistic BEB model for all subshells.

- Start here: `examples/example-Mf.jl`, dated **26-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Mg) Test  computeResonanceStrength(ForDielectronicRecombination(), ...)  for H-like Ne (1s->2l) and Li-like Sc (2s->2p).

- Start here: `examples/example-Mg.jl`, dated **27-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Mh) Test  computeChargeStateDistribution(ForStepwiseDecay(), ...)  for the K-hole decay cascade in Ne and Ar.

- Start here: `examples/example-Mh.jl`, dated **28-Jun-2026**
- Keywords: ForPedestrians, quick start, beginner

### Mi) Test  computeProperties(HyperfineStructure(), ZeemanStructure(), ...)  for HFS A,B and Lande g_J.

- Start here: `examples/example-Mi.jl`, dated **18-Aug-2026**
- Keywords: ForPedestrians, quick start, beginner

---

## N) Examples, tests & development of empirical computations (module Empirical)

### Na) Apply & test Empirical.bindingEnergy (Shell & Subshell), ionizationPotential, and totalEnergy for neutral atoms and ions.

- Module: `Empirical.bindingEnergy`
- Start here: `examples/example-Na.jl`, dated **15-Jul-2026**

### Nb) Apply & test Empirical.photoemissionEinsteinA (ScaledHydrogenic & UsingJAC) and photoexcitation/deexcitation plasma rates.

- Module: `Empirical.photoemissionEinsteinA`
- Start here: `examples/example-Nb.jl`, dated **15-Jul-2026**

### Nc) Apply & test Empirical.photoionizationCrossSection and photorecombinationCrossSection and the corresponding plasma rate coefficients.

- Module: `Empirical.photoionizationCrossSection`
- Start here: `examples/example-Nc.jl`, dated **16-Jul-2026**

### Nd) Apply & test Empirical.impactExcitationCrossSection (Van Regemorter) and impactIonizationCrossSection (Lotz) and their plasma rate coefficients.

- Module: `Empirical.impactExcitationCrossSection`
- Start here: `examples/example-Nd.jl`, dated **17-Jul-2026**

### Ne) Apply & test Empirical.threeBodyRecombinationPlasmaAlpha and the recombination (capture) channels of an ion.

- Module: `Empirical.threeBodyRecombinationPlasmaAlpha`
- Start here: `examples/example-Ne.jl`, dated **16-Jul-2026**

### Nf) Apply & test Empirical.stoppingPower of electrons in a plasma (Bohr1913, Bethe1931, KozmaFranson1992, Axelrod1980).

- Module: `Empirical.stoppingPower`
- Start here: `examples/example-Nf.jl`, dated **17-Jul-2026**

### Ng) Apply & test Empirical.tunnelingIonizationRate, the (quasiclassical) ADK1986 tunneling ionization rate.

- Module: `Empirical.tunnelingIonizationRate`
- Start here: `examples/example-Ng.jl`, dated **17-Jul-2026**

### Nh) Apply & test Empirical.chargeExchangeCrossSection/PlasmaAlpha (OverBarrierModel1980, NiehausScaling1986).

- Module: `Empirical.chargeExchangeCrossSection`
- Start here: `examples/example-Nh.jl`, dated **17-Jul-2026**

### Ni) Apply & test Empirical.neutralizationReducedRate/deExcitationReducedRate, the Belyaev-Yakovleva (2017) model for inelastic A/H,H^- collisions.

- Module: `Empirical.neutralizationReducedRate`
- Start here: `examples/example-Ni.jl`, dated **23-Jul-2026**

### Nj) Apply & test Empirical.dielectronicRecombinationPlasmaAlpha (Arnaud1985DR), for H-like and He-like ions.

- Module: `Empirical.dielectronicRecombinationPlasmaAlpha`
- Start here: `examples/example-Nj.jl`, dated **24-Jul-2026**

### Nk) Apply & test Empirical.excitationAutoionizationPlasmaAlpha (Arnaud1985EA), additive to the Lotz EII rate.

- Module: `Empirical.excitationAutoionizationPlasmaAlpha`
- Start here: `examples/example-Nk.jl`, dated **24-Jul-2026**

### Nl) Apply & test Empirical.generateInelasticHChannels, molecular-symmetry channel correlation for A+H(H^-) collisions.

- Module: `Empirical.generateInelasticHChannels`
- Start here: `examples/example-Nl.jl`, dated **25-Jul-2026**

---

## O) Examples, tests & development of particle-scattering computations (module ParticleScattering)

### Oa) Elastic scattering of electrons at atoms: Dirac partial waves, DCS, Sherman function, transport cross sections.

- Start here: `examples/example-Oa.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Ob) Elastic electron scattering: the high-energy Rutherford limit as an absolute check of the cross sections.

- Start here: `examples/example-Ob.jl`, dated **17-Aug-2026**

### Oc) Elastic electron scattering: the Sherman function and its growth with the nuclear charge.

- Start here: `examples/example-Oc.jl`, dated **17-Aug-2026**

### Od) Elastic electron scattering: the transport cross section of He against measured values.

- Start here: `examples/example-Od.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

### Oe) POSITRON elastic scattering at He beside the electron case.

- Start here: `examples/example-Oe.jl`, dated **17-Aug-2026**

### Of) ONE-PHOTON annihilation of a positron with a bound electron: e^+ + |i(N)> --> |f(N-1)> + photon.

- Start here: `examples/example-Of.jl`
- Date: **unknown/undated** -- NOT VERIFIED; two checks fail, see the report below.

---

## P) Examples, tests & development of photon-scattering computations (module PhotonScattering)

### Pa) BOUND-FREE PAIR CREATION: gamma + |i(N)> --> |f(N+1)> + e^+, the created electron being captured into a bound orbital.

- Start here: `examples/example-Pa.jl`
- Date: **unknown/undated** -- NOT VERIFIED; the gauges disagree, see the report below.

### Pb) RAYLEIGH scattering: gamma + |i> --> |i> + gamma, second order over an explicit gMultiplet of intermediate levels.

- Start here: `examples/example-Pb.jl`, dated **28-Aug-2026**

### Pc) RESONANT INELASTIC scattering (RIXS): the same second order, tuned ONTO an intermediate level with a width.

- Start here: `examples/example-Pc.jl`
- Date: **unknown/undated** -- the Lorentzian IS reproduced to 0.6 %, but the absolute cross sections rest on

### Pd) FORM-FACTOR approximation for Rayleigh and Compton: the cheap limit, and the Thomson absolute check.

- Start here: `examples/example-Pd.jl`, dated **22-Aug-2026**

### Pe) TWISTED (Bessel) beams: the vortex node on axis, and the plane-wave limit, at finite impact parameter.

- Start here: `examples/example-Pe.jl`
- Date: **unknown/undated** -- see the report below.

---

## Q) Muonic atoms

### Qa) Muonic atoms: muon binding energies, muonic X-rays, and the two screenings.

- Start here: `examples/example-Qa.jl`, dated **25-Aug-2026**

---

## S) Stark / line-broadening (semiempirical) computations

### Sa) Apply & test the StarkBroadening module: the electron-impact Stark width of a line, via Dimitrijevic & Konjevic's Simplified Modified Semiempirical (SMSE) method -- needs only the line's own two levels (effective quantum number, active subshell), no perturber search.

- Module: `StarkBroadening`
- Start here: `examples/example-Sa.jl`, dated **21-Oct-2026**

---

## T) Three-body recombination

### Ta) Apply & test the ThreeBodyRecombination module: the rate coefficient alpha^TBR(T_e) via the Saha/Milne detailed-balance relation to electron-impact ionization, from JAC's own computed level energies and a selectable ImpactIonization model (BEB/BED/relativistic/...).

- Module: `ThreeBodyRecombination`
- Start here: `examples/example-Ta.jl`
- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.

---

*Mechanically generated by `tools/buildCapabilityRegister.jl` on 17-Nov-2026 from commit `c64609e1`. Regenerated once per MAJOR/MINOR release (`CLAUDE.md` Rule 17) -- never hand-edited in between.*
