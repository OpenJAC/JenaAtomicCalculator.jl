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

Files outside the JAC project root are **read-only**.
Information exchange between projects is always read-only.

## "Last successful" health contract (Rule 7)

A `# Last successful:  DD-Mon-YYYY` date is written into an example file **only** after the output
has been verified for physical consistency — not merely "it ran without error." Zero rates, wrong
units, or clearly wrong magnitudes mean the date stays blank. A blank date is a regression signal.

Whenever a branch's date is left as `# Last successful:  unknown ...` (i.e. not yet verified), add a
`# Last visit:  DD-Mon-YYYY` line directly above it, updated every time that branch is touched or
run, even without a successful verification. This makes it easy to see later how recently a
not-yet-verified branch was last worked on, independent of whether it ever became verifiable.

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
cd <JAC-root>/test
julia --project=.. runtests.jl
```

Never run `julia --project=. test/runtests.jl` from the JAC root.

## The radial box must be matched to the orbitals (Rule 12)

A box that is far **too large** starves the B-spline basis exactly as badly as one that is too small: the
number of splines is fixed, so a box much wider than the orbitals spends them on empty space and the
high-`n` members of each symmetry are misrepresented -- often by returning a *different state* rather than
an inaccurate one. Symptoms look like angular-momentum or `kappa`-sign bugs and are not.

A hydrogenic orbital `(n,l)` has its outer turning point at `r_+ = (n^2/Z_eff)(1 + sqrt(1 - l(l+1)/n^2))`,
and a box of roughly `2.5 r_+` is a good choice. Use the **effective** charge: an outer electron of a
near-neutral system sees `Z_eff ~ Z - N + 1`, not `Z`.

Between 9 and 10-Aug-2026 this single fact explained four separate "bugs": the E3 1s-5f_7/2 rate that came
out 1000x too small, the Zeeman N1 `kappa <= -3` failure open since July, the `MultipolePolarizibility`
perturbers that were blamed on a B-spline defect, and three successive wrong diagnoses of my own. Before
attributing anything to the angular machinery, **re-run on a box matched to the orbitals**.

Two guards now exist, and neither is sufficient alone:
- `Bsplines.checkGridRepresentation` -- solves the point-nucleus Dirac problem on the given grid and
  compares with `Basics.computeDiracEnergy`. It tests orbitals at the FULL nuclear charge, so it is
  conservative for inner orbitals and **too lenient for diffuse outer ones**.
- `Bsplines.checkOrbitalConsistency` -- run on the FINAL orbitals; the two spin-orbit partners of a subshell
  must agree in mean radius. This is what catches a wrong state. Both honour `AsfSettings.gridStopper`.

## Parking a module (Rule 13)

A module that is deliberately set aside carries a `STATUS, <date>:` block at the end of its module
docstring -- in the repository, not only in memory -- containing four things:
1. the decision and whose it was ("POSTPONED BY THE MAINTAINER"), so a later reader does not treat it as an
   oversight to be helpfully fixed;
2. an explicit instruction not to pick it up without asking;
3. what is concretely known, in enough detail to resume from;
4. a warning against misreading the state -- in particular, which recent edits are *not* progress.

Nothing enforces this; it is documentation. Where a parked module can still be called and would produce
untrustworthy numbers, raise instead, following the `corePolarization.doApply` pattern in
`module-PhotoEmission.jl`: an error that explains, rather than a wrong result or a crash at a random
undefined name.

## This file is public: nothing private leaks out (Rule 14)

`CLAUDE.md` is tracked and published with the repository. Anything written here reaches everyone who
clones JAC, so it must contain **only** what belongs to JAC itself:

- **No absolute paths and no usernames.** Write `<JAC-root>/test`, not a real home directory.
- **No other projects.** Commands, file names or conventions belonging to a different repository, a
  manuscript in preparation, or a collaborator's work do not go here even as a passing mention -- a
  dangling reference to a private paper source is itself a leak.
- **No third-party names**, no unpublished results, no collaborator correspondence.
- **No credentials or machine specifics.**

The working record that *does* hold such things -- findings in progress, wrong turns, colleagues' names,
unpublished plans -- lives outside the repository in the assistant's own memory directory, and stays there.
Keep the two separate: this file is the curated, normative document; that one is the notebook.

## Code hygiene for a source module (Rule 15)

The layout rules a `src/module-*.jl` file should satisfy. They are applied by the `/hygiene` command below;
`/audit` reports on them without changing anything. None of this is about being perfect -- it is about a
reader finding what they expect where they expect it.

**Order within the file**

1. **Structs and their associated methods come first** -- each struct together with its constructors and its
   `Base.show` -- and only then the compute-, determine- and display-methods.
2. **Those methods appear in alphabetical order**, and each is separated from the next by **two blank lines**.

**Docstrings**

3. **Every function and method has its OWN docstring**, even where two of them look nearly alike. A shared
   docstring covering several methods is not acceptable: a reader arriving at one method must find its
   description attached to it.
4. Each docstring opens with the signature -- `` `Module.functionName(argument::Type, ...)` `` -- character
   for character as the `function` line below it, then a line starting with `...` that says what the function
   is **good for**, and it makes the **return value explicit**: "a `value::Type` is returned".
5. **140 characters is the line width to USE, not merely a ceiling not to exceed.** It applies to docstrings,
   comments and code alike. Prose that wraps at 100 or 110 wastes a third of the line and reads worse for it, so
   a paragraph is re-flowed to fill the width; a sentence that fits on one line is never split; and a wrapped
   signature is put back on one line whenever the whole of it fits.

   This applies to the `+ field ::Type   ... description` argument tables too, and they are the main reason for
   the rule: the `+ name ::Type` part is an aligned table and stays exactly as it is, but the description after
   the `...` is ordinary prose that must run to 140 before wrapping, with any continuation indented to the same
   column as the `...`. A one-line description broken over two lines is precisely what this rule exists to stop.

   Only two things keep their own layout, because there the line breaks carry meaning: the keyword list in a
   copy-constructor docstring, and indented example or code blocks.

**Comments**

6. Inline comments start with a single `#`. Never `##`, `###` or deeper.
7. **A `#####...#####` rule that is nothing but hashes is an OPTICAL SEPARATOR and is KEPT.** JAC uses it,
   usually two lines of 129 hashes together, to divide a file into blocks, and that is worth having. What must
   go is a `#####` line that CARRIES TEXT -- `### The display section ####...` -- or a banner block whose
   delimiters enclose comment lines: if such a block says something real it belongs in a docstring, and if it
   only labels the code below it, the code should be findable without it. A detector for this rule must
   therefore test whether the line contains anything besides `#` and whitespace, and must not simply count
   `#####`. Getting this backwards deletes exactly the separators the maintainer wants: on 15-Aug-2026 a
   hygiene pass removed 22 pure separators and zero text-carrying banners.
8. **No empty `#` lines.** A blank line separates blocks of code better than a bare comment marker.
9. Comments are for a reader who knows Julia and atomic physics but not this file's history. **Remove
   historical remarks that can no longer be checked** -- references to code that has since been deleted, to
   an older variant, or to a commit. A docstring says what the code does, not what it once differed from.
   (The `#== ... ==#` form for a large deliberately-disabled block stays; it is a different thing.)

**Bodies**

10. Every compute-method **ends with a clear `return( ... )` line, separated from the body by a blank line**,
    so the result and its type are visible without reading upwards. Very short functions may omit the blank
    line.
11. Top-level `function` lines and their docstrings start at column 1; 4-space indent inside.
12. Existing column alignment in copy-constructors is preserved exactly (see the style section above).

**A hygiene pass changes no behaviour.** That is what makes it verifiable, and the verification is not
optional: the printed output of a real case must be **byte-identical** before and after, and the test suite
must be unchanged with no approved reference re-approved. If a number moves, the pass did something it should
not have.

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
cd <JAC-root>/test && julia --project=.. runtests.jl
```
**Must be run from the `test/` subdirectory** (see Rule 11).
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

### /hygiene ModuleName
**Apply the Rule 15 code-hygiene rules to ONE source module.**
The acting counterpart of `/audit`, which only reports. One module per invocation, as always.

1. Capture a **behaviour baseline** first: run one real case through the module and save its printed output.
   Without this the pass cannot be verified, so it is not optional.
2. Reorder: structs and their methods first, then the remaining methods alphabetically, two blank lines apart.
3. Fix comments (`##` -> `#`, drop TEXT-CARRYING `#####` banners but KEEP the pure hash-only separators, and
   drop empty `#` lines), then **re-flow docstring prose and
   comment blocks to USE the 140-column width** (Rule 15.5) -- filling the line, not just staying under it --
   including the descriptions in the `+ field ::Type ... text` tables, and leaving only copy-constructor keyword
   lists and indented examples alone. A re-flow is verified by
   comparing the **word sequence** before and after: it must be identical, since a re-flow may change only where
   the line breaks fall. That check is what separates re-flowing from rewriting.
4. Give every method its own docstring with signature, purpose and explicit return type; **remove historical
   remarks that can no longer be checked**. Two shapes are worth grepping for, because both are left behind by
   a retirement and both survive as text that no longer says anything: a docstring that describes itself --
   `... as Module.thisSameFunction, but ...`, which was a real contrast only while the other version existed --
   and a paragraph explaining what an earlier form *used to* do.
5. Add the trailing blank line before each `return( ... )` where the function is long enough to warrant it.
6. Verify: the saved output must be reproduced **byte for byte**, and `/test JAC` must be unchanged with no
   approved reference re-approved. Report both.
7. Report what was changed as counts (comments, banners, docstrings added, methods reordered) plus any
   historical remark removed, since that last one is a judgement and the maintainer may disagree.

Steps 2-5 are mechanical and are shown as a representative diff before being applied broadly, per the
collaboration rule on bulk changes.

**Before starting, switch the session to auto-accept** (shift+tab in the CLI and in the VSCode extension).
A hygiene pass is dozens of small scripted steps -- reorder, re-flow, re-measure, re-run the baseline -- and
confirming each one individually costs more attention than reading the final diff does. This is safe here for
a specific reason, not as a general habit: the pass is verified by a byte-identical baseline and an unchanged
test suite, so a wrong step is caught by the verification rather than by the confirmation prompt. The two
prompts still worth keeping are the ones the collaboration rules protect -- the commit and the push.

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

### /priority
**Display the working priority list, and drop from it whatever has been closed.**
The list itself is not kept here — it names collaborators and unpublished work, which Rule 14 keeps out of
this file — but lives in the assistant's memory directory as the single file that the command reads and
rewrites.

1. Read the priority list and print it **in full**, in its stored order, with each item's stored number.
2. For each item, judge whether it is still open. An item is **closed** when the defect is fixed in `src/`,
   the validation has been done, or the maintainer has said it is done — not when it merely looks stale.
3. **Remove every closed item from the stored list, and NEVER RENUMBER.** The numbers are permanent
   identifiers: a closed item's number is retired and the list is left with a gap, which is correct and must
   not be tidied away. Renumbering makes "item 8" mean something different from one day to the next, so that
   neither side can refer to an item reliably. Record each closure in the list's CLOSED section, by number,
   and say which items were removed and why, in one line each.
4. Items carrying a `verify` or `memory only` marker have not been re-checked against `src/`. Do not remove
   one on the strength of the marker alone; either check it or leave it, and say which.
5. Never silently reintroduce an item the maintainer has dropped. The list carries its own dropped-register
   for exactly this.

The command only reads, prints and prunes. It never starts work on an item.

### /application
**Display the application list, and drop from it whatever has been delivered.**
The third of the three working lists, beside `/priority` (defects and open work in `src/`) and the competition
list (capability other codes have). Like the priority list it is not kept here — it names collaborators,
unpublished work and correspondence, which Rule 14 keeps out of this file — but lives in the assistant's memory
directory as the single file the command reads and rewrites.

**What belongs on it.** Physics that is driven from OUTSIDE the code: a colleague is waiting, a paper is in
preparation, a reply is owed, a measurement has appeared that we said we would compare against. The
distinguishing test is that neither of the other two lists' questions decides anything — such an item can be
scientifically minor and still be the most urgent thing here because someone is blocked on it, and it can be
fine physics and wait a year because nobody is. It is therefore ordered by **the strength of the outside
obligation**, not by severity, and the field that matters is *who is waiting, and since when*.

1. Read the application list and print it **in full**, in its stored order, with each item's stored number.
2. For each item, judge whether the obligation is discharged — the run made, the reply sent, the collaborator
   answered — not merely whether the physics still interests us.
3. **Remove every discharged item, and NEVER RENUMBER.** Its number is retired and the list keeps the gap.
4. Its numbers carry an **`A` prefix (A1, A2, …), deliberately not integers**, so that "item 8" always means the
   priority list and "A3" always means this one. The five items it started with arrived from the priority list
   as 36–40, and those integers stay retired there rather than being reused.

The command only reads, prints and prunes. It never starts work on an item.

### /summarize
**End-of-session structured summary.**
Lists: modules changed (and what kind of change), example branches tested (with dates),
open items (blank dates, known bugs, planned next steps).
Important findings are saved to the memory system.

### recall commands
Print the command list above (this section). No steps performed.

## Reference file update protocol

When a code change legitimately alters numerical output (a physics fix or new feature):
1. The new `.sum` output must be verified against an independent source before updating.
2. Updating a `test/approved/*.sum` reference file is a deliberate editorial act — never do it
   automatically or silently.
3. A purely structural change (rename, refactor) must pass all tests unchanged. If a test fails
   after such a change, something went wrong.
