# JAC Coding Conventions for Claude

This file is read by Claude Code at the start of every session. It records the coding conventions
and collaboration rules for the JenaAtomicCalculator (JAC) codebase. Follow these strictly.

## Collaboration rules

- **Never change more than one module (`.jl` file) per task.** Always confirm with the user before
  moving to the next module.
- **All changes must be confirmed by the user before committing.**
- **Never push to the remote without explicit user approval**, even if the user approved a commit.
- **Never make mechanical bulk changes** (regex replacements, automated refactors) without first
  showing a representative sample of the before/after diff and getting approval.
- Commits should group one logical change only. Prefer infrequent, well-described commits.

## Code style — preserve as-is unless told otherwise

### Column alignment in copy-constructors
The keyword copy-constructor pattern in JAC aligns columns deliberately. Do NOT destroy this
alignment when making replacements. Example of the correct style:

```julia
if  Z         == nothing   Zx          = nm.Z           else   Zx          = Z          end
if  Q         == nothing   Qx          = nm.Q           else   Qx          = Q          end
if  Omega     == nothing   Omegax      = nm.Omega       else   Omegax      = Omega      end
```

The `== nothing` test, the assignment targets, the `else`, and the `end` all occupy fixed columns.
A mechanical replacement that shortens or lengthens the test expression will break this layout.

### Naming conventions
- **Type names**: CamelCase (e.g. `Nuclear.Model`, `AbstractPlasmaScheme`, `LineShiftScheme`)
- **Function names**: camelCase or lowercase (e.g. `computeAmplitudes`, `setDefaults`)
- **Variable names**: camelCase, short but descriptive
- **Global constants**: ALL_CAPS_WITH_UNDERSCORES (e.g. `GBL_FRAMEWORK`, `FINE_STRUCTURE_CONSTANT`)
- **Keyword copy-constructor variables**: append `x` to the parameter name (e.g. `Z` → `Zx`)

### Comment style
- `#` — inline explanatory note on a code line
- `##` — section separator or explanatory block comment (used liberally)
- `###` — inline option list (lists the valid alternatives for a setting)
- `#== ... ==#` — large disabled code block (with a reason written on the opening line)
- Never use `##x` — those were dead-code markers and have been removed

### Indentation and spacing
- 4-space indentation throughout
- `if  condition` uses **two spaces** between `if` and the condition (JAC convention)
- Single-line `if ... else ... end` is the standard form for copy-constructor guards
- **Always two blank lines** between the closing `end` of one function/method and the opening
  `"""` docstring of the next. Never collapse to one blank line when editing.
- **Top-level functions and docstrings start at column 1** (no indentation). This applies to all
  `function Module.name(...)` definitions and their preceding `"""..."""` docstrings in `.jl`
  source files. Exception: code inside `include`-d files that is itself inside a module block
  follows the 4-space indent rule like all other nested code (loops, `if-else-end`, etc.).

### Docstring style
Every exported function and struct has a docstring in triple-quoted `"""..."""` form with:
- A header line: `` `Module.FunctionName(args)` ``  followed by a description
- Argument documentation as `+ fieldName ::Type ... description` bullet points
- Examples or cross-references where helpful
- **The function name, argument names, and types in the docstring header must be character-for-character
  identical to those in the `function` line below it.** Never write a supertype or abstract type in the
  header where the function itself uses a concrete type (e.g. write `::Emission`, not `::AbstractEmissionKind`,
  if the method is dispatched on `Emission`).

### Settings structs
Every process/property module has a `Settings` struct with:
- `struct Settings <: AbstractPropertySettings` (properties) or `<: AbstractProcessSettings` (processes)
- A `Settings()` zero-argument default constructor
- A `Settings(set::Module.Settings; kw...)` keyword copy-constructor with aligned columns
- A `Base.show(io::IO, settings::Module.Settings)` method

The keyword copy-constructor guards use `isnothing(field)`, **not** `field == nothing`:
```julia
if  isnothing(calcK)            calcKx          = set.calcK          else   calcKx          = calcK          end
if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end
```
Keyword parameters in the signature use `Union{Nothing,T}=nothing` (e.g. `calcK::Union{Nothing,Bool}=nothing`).

## Directory contract (Rule 6)

| Directory | Purpose | Write? |
|---|---|---|
| `src/` | JAC source modules only — no data, no examples | Yes, one module per task |
| `examples/` | Example scripts `example-Xx.jl` | Yes |
| `demos/` | Pluto notebooks | Yes |
| `work/` | Running JAC; `.sum` output files | Yes |
| Any other project directory | Reference only | **Never** |

Files outside the current project (`/home/fritzsch/fri/JAC.jl/`) are **read-only**.
Information exchange between projects is always read-only.

## "Last successful" health contract (Rule 7)

A `# Last successful:  DD-Mon-YYYY` date is written into an example file **only** after the output
has been verified for physical consistency — not merely "it ran without error." Zero rates, wrong
units, or clearly wrong magnitudes mean the date stays blank. A blank date is a regression signal.

## New ForPedestrians function template (Rule 8)

Every new function added to `module-ForPedestrians.jl` must include, in order:
1. A docstring with a `Simplified call:` block
2. A `Basics.displayConfigurations(stdout, configs, details = "...")` call
3. A hint block (`sa = "\n* Compute ..." * ...`) listing assumptions and optional arguments
4. A `printout=false` suppression path via `redirect_stdout(devnull)`
5. A corresponding entry in `examples/examples.jl` under the M branch
6. A corresponding `example-Mx.jl` file with at least two `if/elseif` scenarios

## Test suite = "Last successful" dates (Rule 9)

The M-branch example files (`example-Ma.jl` through `example-Mf.jl` and beyond) are the active
test suite for `module-ForPedestrians.jl`. Before any structural refactor of a ForPedestrians
function, re-run the corresponding example and confirm the date. A date going blank after a change
is a regression.

## Naming conventions for files and outputs (Rule 10)

- Example files: `example-Xx.jl` where X is the branch letter (A–Z) and x is a–z
- Summary output files: `zzz-<ModuleName>.sum` written to `work/`
- New source modules: `module-<Name>.jl` in `src/`, following existing naming exactly
- Diagnostic/scratch scripts in `work/`: `diag-<topic>.jl`

## Test suite working directory (Rule 11)

`test/runtests.jl` **must always be started from the `test/` subdirectory**, not from
the JAC root.  Every `TestFrames` include file opens summary files with a relative path
(e.g. `"test-Einstein-new.sum"`), so the working directory determines where those files
land.  Running from the root scatters `test-*.sum` and `zzz-*.sum` files across the
JAC root directory.  Running from `test/` places them correctly next to `test/approved/`.

Correct invocation:
```
cd /home/fritzsch/fri/JAC.jl/test
julia --project=.. runtests.jl
```

Never run `julia --project=. test/runtests.jl` from the JAC root.

## Commands

A **command** is a named sequence of steps I execute and then summarize. The leading `/` is optional —
`test Mc` and `/test Mc` are identical. Type `recall commands` at any time to get this list again.

### /test Mx
**Run one example scenario and verify it.**
1. Read `examples/example-Mx.jl`; identify the branch with `if true` (the active scenario).
2. Either the user runs `include("../examples/example-Mx.jl")` in their Julia session and pastes the output,
   or I run it via bash (slower — full JAC compilation).
3. Check the output for physical consistency: correct quantum numbers, right-order magnitudes,
   partial sums matching totals, gauge agreement in expected range.
4. If correct: update `# Last successful:  DD-Mon-YYYY` in the file. If wrong: diagnose before touching code.
5. Report a one-paragraph verdict.

### /test JAC
**Run the general JAC test suite to verify the codebase is not broken.**
Executes `test/runtests.jl` via:
```
cd /home/fritzsch/fri/JAC.jl/test && julia --project=.. runtests.jl
```
**Must be run from the `test/` subdirectory** (see Rule 11 below).
This runs all active `@testset` blocks (methods, structs, evaluations, representations,
amplitudes, properties, processes, cascades). Produces a Pass/Fail summary per testset.
A Fail is reported with the failing test name and error before anything else is done.
Note: some tests are deliberately disabled in `runtests.jl` (marked `##`); those are not failures.
The M-branch `example-Mx.jl` files are **separate** from this test suite.

### /push JAC
**Test, then commit to git.**
1. Run `/test JAC` first. Any failure blocks the push.
2. Collect all modified files (`git status`), draft a commit message, show it to the user.
3. **Wait for explicit user confirmation** before `git add` / `git commit`.
4. After commit: **wait again** for explicit push approval before `git push`.
   (CLAUDE.md rules "never commit / never push without explicit approval" are never bypassed.)

### /audit ModuleName
**Check one source module for convention violations.**
Read `src/module-<Name>.jl` and report anomalies as a bullet list (no changes made):
missing or malformed docstrings, functions not at column 1, misaligned copy-constructor columns,
`##x` dead-code markers, naming violations, missing `Settings()` zero-arg constructor,
`field == nothing` guards instead of `isnothing(field)`.

### /diff
**Summarize recent changes.**
Run `git log` and `git diff` since a reference point (last commit, named date, or "since we started X").
Report: which files changed, what kind (bug fix / new feature / display fix),
and whether any `# Last successful:` dates went blank.

### /diff src
**List changed source files with line counts.**
Runs `git diff --stat HEAD -- src/` and reports only the files under `src/` that differ from the
last commit, together with the number of lines added (+) and removed (−) in each.
Output is a compact table — no diff content shown.
A different baseline can be specified: `diff src <commit-sha>` or `diff src <branch>`.

### /diff src details
**List changed source files and show the actual diffs.**
Runs `git diff HEAD -- src/` and reports the full unified diff for every changed file under `src/`.
Each file section shows the exact lines added (+) and removed (−).
Same baseline rules as `/diff src`.

### /scan src
**Check all source modules for structural anomalies.**
Walk every `src/module-*.jl` file and report: non-source content (data, examples),
modules missing a `Settings` struct, functions not starting at column 1, missing `##` section separators.
Produces a short checklist only — no changes made.

### /summarize
**End-of-session structured summary.**
Lists: modules changed (and what kind of change), example branches tested (with dates),
open items (blank dates, known bugs, planned next steps).
Important findings are saved to the memory system.

### /update MS
Already defined in the article's CLAUDE.md: compile the pedestrian review article PDF via
`pdflatex -shell-escape -interaction=nonstopmode b26.atoms-pedestrian-approach.tex`.

### recall commands
Print the command list above (this section). No steps performed.

## Reference file update protocol

When a code change legitimately alters numerical output (a physics fix or new feature):
1. The new `.sum` output must be verified against an independent source before updating.
2. Updating a `test/approved/*.sum` reference file is a deliberate editorial act — never do it
   automatically or silently.
3. A purely structural change (rename, refactor) must pass all tests unchanged. If a test fails
   after such a change, something went wrong.
