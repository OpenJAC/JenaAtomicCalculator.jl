# How JAC is developed

*Linked from the documentation's front page and from the README.*

## AI assistance, stated plainly

Most changes to JAC are now made with the assistance of an AI coding assistant, working under a written set of
rules kept in [`CLAUDE.md`](https://github.com/OpenJAC/JenaAtomicCalculator.jl/blob/master/CLAUDE.md) in this
repository. Every change was proposed to, and approved by, the maintainer before it was committed, and commits
made this way carry a co-author trailer, so the record is machine-readable and complete; the table below gives
the date of the earliest such commit rather than a hand-remembered one:

<!-- BEGIN generated: development record -->
| | |
|---|---|
| commits in the whole history | 1452 |
| of those, carrying an AI co-author trailer | **1119** |

| year | commits | carrying the trailer |
|---|---|---|
| 2019 | 102 | 0 |
| 2020 | 56 | 0 |
| 2021 | 22 | 0 |
| 2022 | 22 | 0 |
| 2023 | 17 | 0 |
| 2024 | 14 | 0 |
| 2025 | 40 | 0 |
| 2026 | 1179 | 1119 |

Of the commits carrying the trailer, the earliest carries author date **9 June 2026**.

*Measured 16 November 2026 at commit `d4cc7939`; regenerate with `julia --project=. tools/measureDevelopmentRecord.jl --write`.*
<!-- END generated: development record -->

The physics judgements, the priorities, and the decisions not to proceed are the maintainer's.

## Why the record is kept this way

The documented risk in AI-assisted scientific computing is not code that crashes but **code that runs while
encoding the wrong physics** — a fitted coefficient that passes every test it was shown, a tolerance quietly
loosened until a wrong answer passes, a verification claimed but not performed. The common remedy in the
literature is the same in every study: a reference the acting party cannot edit.

JAC's rules were arrived at independently and say, in part:

- **One source module per change**, and no mechanical bulk edit without an approved sample diff.
- **A pass that claims to change no behaviour must reproduce a real case's printed output byte for byte**, and
  no approved test reference may be re-approved in order to make the suite pass.
- **Three date markers on an example branch, of which only two may be written by a tool.** That a branch ran
  (`Last visit`) and that its own quoted numbers came back (`Last reproduced`) are mechanical facts. That the
  *physics* was judged correct (`Last successful`) is recorded only by a person. The reason is a case in this
  repository: a branch whose every number reproduced exactly had its date withdrawn, because the comparison
  being made was between two different physical quantities. **Reproduction proves stability, not correctness.**

## What this means for a user of JAC's numbers

It means the same thing it would mean for any atomic-structure code: a computed number is worth what its
verification is worth. Where an example branch carries a `Last successful` date, a person has judged the
physics on that date, and the branch quotes the numbers that were judged. Where it does not, the branch may
still run perfectly well and has not been through that judgement. The distinction is deliberately visible, and
it is the most useful thing in this documentation.

Independent comparison against other codes and against measurement is the check that neither a test suite nor
a development rule can provide, and it is what the published applications of JAC exist to do.

## Keeping this page honest

The table above is not hand-typed. It is re-derived from `git log` on the repository itself by
[`tools/measureDevelopmentRecord.jl`](https://github.com/OpenJAC/JenaAtomicCalculator.jl/blob/master/tools/measureDevelopmentRecord.jl),
the same way `docs/src/scf-routes.md`'s accuracy table and `CAPABILITIES.md` are re-derived rather than
maintained by hand — run once per release (`CLAUDE.md`, Rule 17) and never hand-edited in between:

```
julia --project=. tools/measureDevelopmentRecord.jl --write
```

A commit counts as AI-assisted if its body carries a `Co-Authored-By:` trailer naming Claude or Anthropic. This
undercounts a commit whose assistant forgot the trailer, which has happened; the number above is therefore a
lower bound, not an exact count, and is presented as one.
