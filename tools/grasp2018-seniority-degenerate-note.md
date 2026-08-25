# A reproducible discrepancy in GRASP2018 two-particle angular coefficients

**Status:** reproducible, attributed with argument but **not** diagnosed at the level of GRASP's source.
Written 25-Aug-2026 from the comparison harness in `tools/diag-grasp-angular.jl`.

This note exists so that anyone maintaining GRASP2018 has a concrete case with numbers, rather than a report that
"another code disagrees". Everything below can be re-run; nothing rests on trusting our side.

## What is affected

Two-electron (Coulomb) angular coefficients for a **diagonal** matrix element of a CSF that shares its occupation,
total `J` and parity with another CSF of the same list — that is, two CSFs differing **only in their intermediate
coupling / seniority**. Coefficients of CSFs uniquely determined by (occupation, `J`, parity) agree with us exactly.

The condition is the one `librang90`'s own ReadMe names: *"configurations whose difference reduced only to the value
of the seniority quantum number"*.

## Reproduction

Generate the CSF list with `rcsfgenerate` from the two configurations

```
1s(2,i)2s(2,i)2p(1,i)
1s(2,i)2s(1,i)2p(2,i)
```

with active set `2s,2p` and `2*J` from 1 to 5, no excitations. This gives 10 CSFs over the peel subshells
`1s, 2s, 2p-, 2p`. Call `RKCO_GG(IC, IR, CORD, 1, 1)` for each pair and read `BUFFER_C`.

Among the CSFs there are two with occupation `1s(2) 2s(1) 2p-(1) 2p(1)`, `2J = 3`, even parity. They are GRASP's
CSFs **6** and **7** in the order `rcsfgenerate` produces. On the **diagonal** of CSF 6, two coefficients differ from
ours; the diagonal of CSF 7 agrees with us on all six of its shared coefficients.

| rank | Slater integral | our value | GRASP2018 |
|---|---|---|---|
| k = 1 | `R^1(2s, 2s, 2p_3/2, 2p_3/2)` | −0.16666667 | −0.04244067 |
| k = 2 | `R^2(2p_1/2, 2p_1/2, 2p_3/2, 2p_3/2)` | −0.10000000 | −0.14472136 |

Both are **exchange** terms. Values are the coefficient of the Slater integral including the
`(-1)^k <a||C^k||c> <b||C^k||d>` factor, i.e. directly comparable with what `SPEAK` deposits.

## Why we attribute it to GRASP rather than to ourselves

Three independent reasons, in increasing weight:

1. **Two independent implementations on our side agree.** JAC's production module and a second module written from
   scratch — different recoupling, different shell tensors, derived rather than transcribed — give identical values.

2. **The condition matches the one GRASP's own documentation flags** as the case a 2009 correction addressed. That
   correction guards only the core section, after `IF (INCOR .LT. 1) RETURN`; the terms above are produced on the main
   path, which the guard does not reach. **This is a pointer for where to look, not a claim that we have read the
   defect.**

3. **It cannot be a basis difference, and this is the decisive one.** The natural objection is that two CSFs differing
   only in coupling span a two-dimensional space in which the two codes might simply carry different bases, making a
   coefficient-by-coefficient comparison meaningless. That objection is testable and fails: the **partner** state
   (GRASP's CSF 7) agrees with ours on **all six** of its shared diagonal coefficients. If the two codes carried
   different bases of one two-dimensional subspace, one basis vector coinciding exactly would force the other to
   coincide too. It does not.

   Equivalently: a trace over the degenerate pair is basis-independent, and ours and GRASP's traces differ — on
   exactly the two orbits above and no others.

## What we have NOT established

- We have not identified the responsible lines in `librang90`. The attribution rests on the three arguments above.
- We have tested one configuration pair. The defect may be broader or narrower than this case.
- We make no claim about the size of any physical consequence. These are angular coefficients; whether a given
  calculation is affected depends on which Slater integrals multiply them and how large those are.

## How to re-run the comparison

`tools/diag-grasp-angular.jl` builds GRASP2018 from a read-only source tree into a scratch directory, drives
`rcsfgenerate`, dumps the coefficients through two small Fortran drivers (`tools/grasp-drivers/`), and compares them
with JAC's, matching CSFs by signature and — where that is degenerate — by their diagonal two-particle coefficients
scored on shared keys. The traps that cost time to find once (library build order, `-std=legacy
-fallow-argument-mismatch`, the `DLAMCH` shim, the `2*J` parity that `rcsfgenerate` requires) are handled there and
documented at the top of the file.
