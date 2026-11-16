# Instructions for coding agents

The instruction file for this repository is **[`CLAUDE.md`](CLAUDE.md)**. Read it before proposing any change;
it is not a style guide but the project's working constitution, and it is binding on assistants and on the
maintainer alike.

This file exists only because `AGENTS.md` is the name several tools look for. It adds no rules of its own.

Three points that are easy to miss and that matter most:

- **One source module per change**, and no mechanical bulk edit without an approved sample diff.
- **`test/approved/*.sum` reference files are never updated to make a test pass.** If a number moves, that is a
  physics question, not a maintenance task.
- **The `Last successful` date on an example branch may only be written by a person.** Reproducing a number
  proves stability, not correctness.

If you are assisting an outside contributor, [`CONTRIBUTING.md`](CONTRIBUTING.md) is the shorter document and
is the one addressed to you.
