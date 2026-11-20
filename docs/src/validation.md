# Validation against another code

Almost every check a code applies to itself is self-referential. JAC's test suite compares against
reference files that JAC itself produced; the cancellation factor and the gauge ratio that flag a
doubtful transition are formed from the very amplitudes they judge; and the `Last successful` date on
an example records one person's reading of JAC's own output. None of that can detect an error that
JAC makes consistently.

This page is the exception: a standing comparison against values that do **not** come from JAC.
It is regenerated mechanically, so nothing on it is typed in by hand.

## How a row must be read

The rows are grouped by **tier**, because the tier decides what a difference *means*. A relative
difference of ``10^{-7}`` is a defect in one tier and entirely expected in another, and a table that
mixes them says nothing at all.

| tier | the reference is | a difference means |
|---|---|---|
| 0 | a closed-form expression — the Dirac energy of a one-electron ion in a point nucleus | the discretization of one of the two sides; it bounds how sharply anything else can be compared |
| 1 | another code solving **the same equations** with **the same configuration state functions** | a defect in one of the two codes |
| 2 | another code using a **different method** | informative, never decisive: the two calculations are not the same calculation |

Within a tier the rows are ordered by the quantity compared, since that is the question a reader
arrives with — *can I trust JAC's transition rates?* — rather than by which partner supplied the
reference. The partner is therefore a column: where two codes cover one quantity, their rows belong
side by side.

## What is pinned, and why

A comparison is only as sharp as the agreement between the two problems being solved, so everything
that could differ is fixed rather than left to each code's defaults:

  * a **point nucleus**, because two finite-nucleus models are never identical and the difference
    would be read as a difference between the codes;
  * the **Coulomb interaction only**, with no Breit term and no QED correction — these are the terms
    a closed-form Dirac energy does not contain, and the terms whose implementations may legitimately
    differ;
  * a **single reference configuration**, so that the list of configuration state functions is one
    function in both codes and no correlation model stands between them. This is what makes a tier-1
    difference a defect, and it is possible only with a partner that uses the same formalism.

The **radial grids are deliberately not made identical**, and that is not an oversight. A
finite-difference integrator carries a value at every radius, while JAC expands each orbital in a
fixed set of B-splines; the same four grid parameters therefore do not describe the same
discretization for the two codes, and a grid well suited to one can starve the other. Each code is
run on a grid adequate for itself, and the comparison rests on the Hamiltonian and the function list
being identical.

## The comparison

<!-- BEGIN generated: cross-code validation -->

#### Tier 0 — against a closed form (this tier measures the floor; it does not grade)

| system | level | quantity | JAC | partner | partner value | rel. diff | verdict |
|---|---|---|---|---|---|---|---|
| H Z = 1 | 1s_1/2 | level-energy | -0.5000066566 | exact-dirac-formula `closed form` | -0.5000066566 | 5.62e-13 | REFERENCE-LIMITED |
| H Z = 1 | 2p_1/2 | level-energy | -0.1250020802 | exact-dirac-formula `closed form` | -0.1250020802 | 5.28e-12 | MEASURED |
| H Z = 1 | 2p_3/2 | level-energy | -0.125000416 | exact-dirac-formula `closed form` | -0.125000416 | 5.28e-12 | REFERENCE-LIMITED |
| Ar17+ Z = 18 | 1s_1/2 | level-energy | -162.7048581 | exact-dirac-formula `closed form` | -162.7048581 | 2.53e-12 | MEASURED |
| Ar17+ Z = 18 | 2s_1/2 | level-energy | -40.72036384 | exact-dirac-formula `closed form` | -40.72036384 | 4.95e-12 | MEASURED |
| Ar17+ Z = 18 | 2p_1/2 | level-energy | -40.72036384 | exact-dirac-formula `closed form` | -40.72036384 | 3.86e-12 | MEASURED |
| Ar17+ Z = 18 | 2p_3/2 | level-energy | -40.54376721 | exact-dirac-formula `closed form` | -40.54376721 | 4.25e-12 | MEASURED |
| Xe53+ Z = 54 | 1s_1/2 | level-energy | -1519.473343 | exact-dirac-formula `closed form` | -1519.473343 | 1.24e-10 | MEASURED |
| Xe53+ Z = 54 | 2p_1/2 | level-energy | -383.790162 | exact-dirac-formula `closed form` | -383.790162 | 9.77e-12 | MEASURED |
| Xe53+ Z = 54 | 2p_3/2 | level-energy | -368.1078699 | exact-dirac-formula `closed form` | -368.1078699 | 3.26e-12 | MEASURED |
| U91+ Z = 92 | 1s_1/2 | level-energy | -4861.197776 | exact-dirac-formula `closed form` | -4861.197904 | 2.62e-08 | MEASURED |
| U91+ Z = 92 | 2p_1/2 | level-energy | -1257.395844 | exact-dirac-formula `closed form` | -1257.395852 | 6.51e-09 | MEASURED |
| U91+ Z = 92 | 2p_3/2 | level-energy | -1089.611416 | exact-dirac-formula `closed form` | -1089.611416 | 1.53e-12 | MEASURED |

#### Tier 1 — same formalism, same CSF list: a difference is a defect

| system | level | quantity | JAC | partner | partner value | rel. diff | verdict |
|---|---|---|---|---|---|---|---|
| He-like Z = 18 | 2J = 0,0 | level-energy | -314.2001643 | GRASP2018 `9006157` | -314.2001644 | 1.82e-10 | AGREES |
| He-like Z = 54 | 2J = 0,0 | level-energy | -3003.180704 | GRASP2018 `9006157` | -3003.180704 | 1.03e-10 | AGREES |
| Be-like Z = 18 | 2J = 0,0 | level-energy | -379.1988155 | GRASP2018 `9006157` | -379.1988162 | 1.65e-09 | AGREES |
| Be-like Z = 54 | 2J = 0,0 | level-energy | -3716.721733 | GRASP2018 `9006157` | -3716.721708 | 6.55e-09 | AGREES |
| Ne-like Z = 26 | 2J = 0,0 | level-energy | -1148.430869 | GRASP2018 `9006157` | -1148.430868 | 6.27e-10 | AGREES |
| Ne-like Z = 54 | 2J = 0,0 | level-energy | -5566.843185 | GRASP2018 `9006157` | -5566.843173 | 2.09e-09 | AGREES |

19 comparisons, 0 marked `DIFFERS`.  Measured 20 November 2026 at commit `b180ad63`; regenerate with `julia --project=. tools/crossCheckBenchmark.jl --write`.
<!-- END generated: cross-code validation -->



The same rows are written to
[`validation-rows.csv`](validation-rows.csv), one line per comparison, carrying each row's model, its
grid, its tolerance and the partner's pinned commit under field names that are expected to stay
stable. That file is there so that another project's harness can read these comparisons, add a column
of its own and check them, without needing anything from us.

## What this does not establish

Every row here tests whether JAC **implements** a given set of equations correctly. None of them can
say whether the equations were the right ones for a particular atom: the choice of correlation model,
of mean field and of basis remains a physical judgement, and agreement with another code is no
substitute for it. A clean table is necessary, and it is a long way from sufficient.
