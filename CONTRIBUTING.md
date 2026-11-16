# Contributing to JAC

Thank you for your interest. JAC is a relativistic atomic-structure code maintained by a very small group, and
that single fact shapes everything below: the scarce resource here is not ideas or code but **the judgement of
whether a computed number is physically right**. These guidelines exist to spend that judgement well.

Please read the two short sections before opening a pull request. They are short on purpose.

## What makes a contribution reviewable

Four things, and a submission that has them can be reviewed in a sitting rather than over months.

1. **One module where possible.** JAC's own development rule is one source module per change, because a reviewer
   has to hold the module in mind while reading it. A submission touching eight modules is eight reviews, and it
   will be slower than eight separate submissions. If the work genuinely cannot be split, say so and say why.

2. **Say which approved references change.** `test/approved/*.sum` files are the reference output of the test
   suite. If your change moves a number in one of them, **do not update the reference** — leave it as it is and
   list the affected files in the pull request. Updating an approved reference is a deliberate editorial act and
   requires a physics decision, which is the maintainer's. A submission that silently re-approves a reference
   cannot be accepted, however good the code is.

3. **Quote the numbers you claim, and say what you checked them against.** "The tests pass" is not a
   verification standard for a physics code; a test suite checks that the code still does what it did. What
   helps a reviewer is a number, a reference value, and the source of that reference — a published table,
   another code, an analytic limit, a sum rule, a gauge comparison. If agreement is partial, say where it fails;
   a stated discrepancy is far more useful than a silent one.

4. **Disclose AI assistance.** See the next section. It is one line and nothing turns on the answer.

## AI assistance

**JAC itself is developed with substantial AI assistance, and says so.** Most changes since January 2026 were
made with an AI coding assistant working under the written rules in [`CLAUDE.md`](CLAUDE.md), with every change
proposed to and approved by the maintainer before commit, and such commits carry a co-author trailer. So this is
not a project that disapproves of the practice.

What we ask of a contributor is therefore not permission but **disclosure and responsibility**:

* **Tell us where an assistant contributed** — tests, derivation, bug-hunting, code, documentation. The
  pull-request template has a field for it. There are no restrictions attached to the answer; it tells a reviewer
  where to read closely, which is useful in exactly the way a note saying "this part is machine-translated from
  Fortran" would be useful.
* **You must be able to explain the code you submit**, including how it relates to the surrounding code. This is
  the one firm requirement, and it is the same one NumPy, PySCF and psi4 state. If you cannot explain a part of
  your own submission, it is not ready.
* **A person submits and answers questions.** A pull request opened autonomously by an agent will be closed.

## The physics date, and what a contribution may claim

Example branches in `examples/` carry three date markers, and the distinction matters for what your pull request
should claim:

| marker | asserts | who may write it |
|---|---|---|
| `Last visit` | the branch was run | anyone, including a tool |
| `Last reproduced` | it ran and its own quoted numbers came back | anyone, including a tool |
| `Last successful` | a **human** judged the physics | a person only, never a tool |

A contribution may legitimately claim the first two. Please do not write a `Last successful` date; that judgement
is made on this side, and leaving it to us is not a formality — reproduction proves stability, not correctness,
and this project has a recorded case of a branch whose every number reproduced exactly while the comparison it
made was the wrong one.

## How a pull request is handled

Three steps, in order, so you know where yours is:

1. **Mechanical** (fast): does it rebase onto current `master`, does the suite run, which approved references
   change, is the AI field filled in?
2. **Scope** (fast): can it be split into independently reviewable pieces? *Being asked to split a submission is
   not a rejection.* Diagnostics and scripts under `tools/` can usually be taken straight away, since that
   directory is neither published nor load-bearing.
3. **Physics** (slow): for every changed number, is the new value right? This is the step that takes real time,
   and it is the reason the first two exist.

Running the test suite, for step 1, must be done from the `test/` subdirectory:

```
cd test
julia --project=.. runtests.jl
```

## Smaller things that are always welcome

An example branch that runs and is dated; a diagnostic under `tools/`; a documentation fix; a clear bug report
with the input that triggers it. A reproducible bug report is worth a great deal here and costs you less than a
patch.

## Questions

Open an issue. A question about whether a change would be wanted is cheaper for both sides than a finished
submission that turns out not to be.
