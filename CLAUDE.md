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
- **Two sessions often run in parallel in this one directory.** The discipline that keeps them from
  blocking or overwriting each other is under `/preparetravel`, "Starting a SECOND session"; read it
  before committing, testing or editing the shared working lists.

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
| `tools/` | Assistant-maintained diagnostics and cross-code comparison harnesses | Yes |
| Any other project directory | Reference only | **Never** |

`work/` is **gitignored**: it is scratch, and anything left there is lost on the next clean. `tools/` is
**tracked** and is the assistant's own -- diagnostics, harnesses, and the small Fortran drivers that let JAC's
results be compared against another code. The maintainer is not expected to read or maintain it; its point is
that such a harness survives a move between machines, which a `work/` file does not.

Files outside the JAC project root are **read-only**.
Information exchange between projects is always read-only.

## One application at a time (Rule 6a)

`apps/` holds the work done FOR PEOPLE -- one subdirectory per collaboration, and inside them job scripts, results,
correspondence, contact details and drafts. It is gitignored and it is not the assistant's to browse. The
subdirectories are separate because the people are.

**A session works in the ONE subdirectory it has been told to work in.** It reads, runs and writes there. It does not
open another `apps/` subdirectory -- not to look for a similar calculation, not for a script that would save time, not
to check how something was done before -- and it carries nothing between them: no numbers, no text, no file, and above
all no name. If material from another application would genuinely help, the assistant ASKS FIRST and says exactly
which directory and what it wants from it.

The point is not tidiness. Each directory belongs to a different collaborator, and information that moves between
them without being asked for is a small breach of their confidence rather than a filing error.

The maintainer names the directory in ordinary words -- "today we work in apps/<name>" -- and nothing else is needed.
The session should still be started on the JAC root rather than on the subdirectory, since application work usually
touches `src/` as well.


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

**A box that is merely TOO SMALL fails in the opposite and more dangerous way: it corrupts progressively, so the
artifact impersonates a physics conclusion.** On 25-Aug-2026 a second-order convergence study of two-photon ionization
gave gauge ratios marching monotonically AWAY from unity as bound intermediates were added -- 0.775, 0.553, 0.396,
0.198 in a 60 a.u. box -- which reads as "the bound intermediate sum does not converge", and was half written up as
exactly that. The same study in a 90 a.u. box converges instead: 0.793, 0.803, 0.873, 0.908. The corrupted continuum
orbital grew worse as the sum grew, so the artifact and the physics carried the SAME signature, and the DIRECTION of
the march is what made it convincing -- a single bad ratio invites suspicion, four in a monotone sequence look like a
trend with a cause.

**Note which failure is loud.** Near threshold `Continuum` REFUSES outright with "enlarge box-size", so the unusable
case announces itself, while the merely-too-small case runs clean and returns numbers. A too-large box produces a
wrong NUMBER; a too-small one produces a wrong CONCLUSION.

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

**AND THE SAME APPLIES TO EVERY COMMIT MESSAGE, WHICH IS THE HARDER HALF.** This file is edited deliberately
and re-read often; a commit message is written once, in the middle of the work, and is then permanent, public and
effectively unfixable -- removing one means rewriting pushed history that other clones already hold. On
31-Aug-2026 a scan of the last 100 commits found exactly one leak, and it had the shape this rule must catch:
the physics was described correctly and a collaborator's FULL NAME was put in brackets beside it, as ordinary
context, by a session that had been working in that person's `apps/` folder all morning.

**So an application is named by its NUMBER and its PHYSICS, never by a person.** Write

    found while starting application A17 (Cf^15+/Cf^17+ q-factors)

and never

    found while starting application A17 (Cf^15+/Cf^17+ q-factors for <a person's name>)

The number is the link to the working record, which is where the name is allowed to live. This holds for the
commit subject and the body alike, for a `.jl` comment or docstring in `src/`, and for a file name -- anything
that git carries. A first name alone is no better than a full one; if a reader could identify the person, it does
not go in.

**The check is one command and takes seconds**, so it is worth running before a release and after any stretch of
work driven from an `apps/` folder: build the token list from the `apps/` directory NAMES, then

    git log -100 --format='%an%n%s%n%b' | grep -inwF -f <that token list>

Read the hits rather than counting them -- `alpha`, `breit`, `plus`, `high` and `thomas` all matched innocently
on 31-Aug (Thomas-Reiche-Kuhn is a sum rule), and the one real hit sat among them. Scan the DIFFS too, not only
the messages: a name committed into `src/` is the worse leak of the two.

The working record that *does* hold such things -- findings in progress, wrong turns, colleagues' names,
unpublished plans -- lives outside the repository in the assistant's own memory directory, and stays there.
Keep the two separate: this file is the curated, normative document; that one is the notebook.

## Code hygiene for a source module (Rule 15)

The layout rules a `src/module-*.jl` file should satisfy. They are applied by the `/hygiene` command below;
`/audit` reports on them without changing anything. None of this is about being perfect -- it is about a
reader finding what they expect where they expect it.

**Order within the file**

0. **`RacahAlgebra-inc-*` IS OUT OF SCOPE FOR HYGIENE ENTIRELY** -- the maintainer's decision, 28-Aug-2026. Do
   not run `/hygiene` on `module-RacahAlgebra-inc-sumrules.jl` or `module-RacahAlgebra-inc-special.jl`, and do
   not re-propose them. Their `#` lines are structural to how the sum rules are written out, not clutter.

1. **Structs and their associated methods come first** -- each struct together with its constructors and its
   `Base.show` -- and only then the compute-, determine- and display-methods.
2. **Alphabetical order applies WITHIN a `#####` banner section, not across the whole file** -- the maintainer's
   instruction, 28-Aug-2026. A banner such as `### RAS: Restricted-Active-Space expansions ####...` marks a block
   whose structs and methods belong to ONE topic, and that grouping carries real information: sorting the file
   end-to-end would scatter a topic's methods among unrelated ones and destroy it. So the banners fix the order of
   the SECTIONS, and alphabetical order arranges the methods INSIDE each section. In a file with no banners,
   alphabetical order applies to the file as a whole, as before.
   Each method is separated from the next by **two blank lines**.

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

6. **A comment line that STANDS ON ITS OWN starts with a single `#`.** Never `##`, `###` or deeper.
   **A comment APPENDED TO A CODE LINE may keep its `##`** -- amended 28-Aug-2026 on the maintainer's
   instruction. The doubled hash is doing work there: it marks the text as an aside to the statement it sits on,
   and converting it buys nothing.
   **What to do with an appended comment instead is to READ IT.** Ask whether it is instructive -- whether it
   tells a reader something the code does not already say. If it does not, DELETE IT rather than reformat it.
   `n = n + 1    ## increment n` is noise whichever hash it uses.
7. **A `#####...#####` rule that is nothing but hashes is an OPTICAL SEPARATOR and is KEPT.** JAC uses it,
   usually two lines of 129 hashes together, to divide a file into blocks, and that is worth having. What must
   go is a `#####` line that CARRIES TEXT -- `### The display section ####...` -- or a banner block whose
   delimiters enclose comment lines: if such a block says something real it belongs in a docstring, and if it
   only labels the code below it, the code should be findable without it. A detector for this rule must
   therefore test whether the line contains anything besides `#` and whitespace, and must not simply count
   `#####`. Getting this backwards deletes exactly the separators the maintainer wants: on 15-Aug-2026 a
   hygiene pass removed 22 pure separators and zero text-carrying banners.
   **AMENDED 28-Aug-2026: a SHORT, INSTRUCTIVE label inside a separator MAY STAY.** The maintainer's own
   example is `####### Radiative Recombination ####...`, and such a line earns its place -- it tells a reader
   scanning the file which block they have reached, which is exactly what these rules are for. The test is
   therefore NOT "does the line carry text" but "does the text help someone find their way". A one- or
   two-word section name inside a separator: KEEP. A sentence, a paragraph, or a delimiter block enclosing
   comment lines: that belongs in a docstring, and goes.
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

**READ THESE RULES WITH TOLERANCE** -- the maintainer's instruction, 28-Aug-2026. The goal is to make the code
**more easily accessible to a potential user**, not to reach a state a detector calls clean. Where a rule and that
goal disagree, the goal wins, and the judgement is worth stating in the report rather than resolved silently. A
mechanical sweep that satisfies every clause while making a file harder to read has failed at the only thing it
was for.

**A hygiene pass changes no behaviour.** That is what makes it verifiable, and the verification is not
optional: the printed output of a real case must be **byte-identical** before and after, and the test suite
must be unchanged with no approved reference re-approved. If a number moves, the pass did something it should
not have.

## What gets published in the documentation (Rule 16)

Agreed 18-Aug-2026, after measuring several candidate rules against the actual code. The problem it solves: the
`@autodocs` blocks in `docs/src/api-*.md` selected content with hand-maintained `Pages` allowlists, so a source
file nobody added was simply absent — and `warnonly = Documenter.except()` in `docs/make.jl` downgraded every
Documenter check to a warning, so the build could not complain. Whole modules were missing from the published
API and nothing said so.

**1. `TestFrames` is never published.** It is the test suite.

**2. Infrastructure modules are published IN FULL**, minus a small name-pattern exclusion for internals:
`*Kernel`, `*Reference`, `*_coefficients`, `*_densities`, `*RouteOf`, plus the primitives `bracket`, `oneJ`,
`oneM`, `j_values`, `m_values`. Measured: this removes **12 functions of 653, 1.8 %**.

   Two rules were tried and REJECTED, and are recorded so they are not retried. *Exported names only* fails
   because JAC's convention is the qualified call: `AngularMomentum` exports 1 name of 70, `InteractionStrength`
   1 of 88, so almost everything would vanish. *Used in an example* fails at function level for the same reason
   — users call `perform(...)`, not internals, so `AngularMomentum` scores 0 of 32.

   The asymmetry decides it: a stray helper in the docs costs a reader one line, a missing main function costs
   them the feature.

**3. A physics module is published if it is USED IN THE CODE of an example file carrying at least TWO
`Last successful` dates.** "Used in the code" means a real `Module.` call with comments stripped — a mention in
a comment does not count, and neither does an `examples.jl` description, which under-counted by more than half
when tested (`DielectronicRecombination` scored 0 by description and 16 by usage). A module with no dated
example is not published.

**4. The `-inc-` split is an authoring convention, not a user-facing boundary**, so it is NOT used for
inclusion: once a module qualifies, ALL of its files are documented and the `Pages` allowlists disappear —
which is what removes the silent-loss bug. The one exception is a file whose own types appear ONLY in UNDATED
example branches: that file has never been used for successful work and is excluded. Measured at branch
granularity (244 dated branches against 168 undated), this excludes **two** files today.

   Applying rule 3 at `-inc-` granularity was tried and REJECTED: only 6 of 67 files would have qualified,
   because examples name the scheme (`Cascade.PhotoAbsorptionScheme`) while the work sits in files whose own
   types are never mentioned.

**The filter is self-correcting, and that is the point.** Nothing is excluded by a hand-maintained list. The day
someone verifies a branch and writes the date, the file publishes itself with no change to `docs/`. Equally, a
module that loses its dated examples drops out. **Re-derive the sets rather than editing them by hand**; a
stale allowlist is the defect this rule exists to remove.

## Releasing a new version (Rule 17)

Written 18-Aug-2026 while doing 0.5.0. It lives here rather than in `docs/NewDocumentationRelease.txt`, which is
gitignored and therefore reaches nobody.

1. **A release is a SNAPSHOT.** Run `git status` and see who else is mid-edit. Work left uncommitted is not in
   the release, and that must be a decision rather than an accident — if a topic is in flight and cannot wait,
   say so in `docs/src/news.md`. A caveat in the release notes is honest; silence is not.
2. **`julia --project=. docs/checkCoverage.jl`** must print OK. It re-derives what Rule 16 says to publish and
   what the API pages actually reference, reading no stored list. If a module qualifies but is undocumented, ADD
   it; do not edit the rule to silence the message.
3. **`julia --project=docs docs/make.jl`** must exit 0. Two failures are worth expecting: a `@example` block
   breaking through ordinary API drift — the blocks are CHAINED, so one break cascades and Documenter stops at
   the first — and `HTMLSizeThresholdError`, where a page passed 1 MB and should be split rather than have the
   threshold raised. Deployment cannot happen locally; Documenter detects it is not in CI and skips.
4. **`docs/src/news.md`**, newest first, selecting what a USER would notice over what the repository noticed —
   0.5.0 drew eight entries from about 220 commits.
5. **Bump `version` in `Project.toml`**, then run the suite from `test/` (Rule 11). No approved reference may be
   re-approved to make it pass.
6. **Commit as `Release X.Y.Z: <one line>`** and push. `documentation.yml` then rebuilds and deploys the site
   automatically, so the docs go live at push, before the registry step.
7. **Register — the maintainer's step.** Comment on the release commit:

       @JuliaRegistrator register

   **WITH RELEASE NOTES IN THE SAME COMMENT.** While the version is 0.x, every MINOR bump is BREAKING by Julia's
   rules, and AutoMerge refuses a breaking release whose notes do not explain the break. This caught 0.4.0 and
   then 0.5.0 again. Derive the list mechanically —

       git log v<PREVIOUS>..HEAD --oneline | grep -iE "retire|un-export|remove|rename|no longer"

   — then keep what a CALLER would notice and drop the internal churn.

   **Re-posting the trigger on an already-registered commit does NOTHING**: on 18-Aug the bot answered the first
   request in 14 seconds and ignored the second for 44 minutes, because version and commit were unchanged. To add
   notes afterwards, register a NEW commit. Do not create the tag by hand — `TagBot.yml` does it once the General
   registry merges, and doing both yields two tags for one release.
8. **Afterwards**, read WHICH CI job failed before reacting: JAC's CI has historically failed on the coverage
   step rather than on the tests.

## Where the sqrt(2j+1) sits: one Wigner-Eckart convention for the whole code (Rule 18)

Agreed 25-Aug-2026, after an inventory found **eighteen** places compensating for this factor across **nine** modules,
in two opposite directions, plus one more where the compensation hides inside an integral.

**First, what is NOT a convention and cannot be normalized away.** A rank-0 coefficient multiplies an ORDINARY matrix
element `<a| o |b>`; a rank-k coefficient multiplies a REDUCED one `<a|| o^(k) ||b>`. Those are different objects. No
choice of factors makes them the same, and any code that treats them alike is wrong rather than unconventional.

**Second, what IS a convention.** `sqrt(2j_a+1)` is exactly the Wigner-Eckart factor between those two objects, so it
may be carried either by the coefficient or by the matrix element. **JAC's choice is GRASP's: the factor sits INSIDE
the coefficient, at EVERY rank.** Three reasons, in order: the rank-0 diagonal coefficient is then literally the
occupation number, so a wrong convention is visible by eye instead of invisible; it matches the published tables and
twenty years of GRASP results; and it is the smaller migration from where the code stands.

**Third, and this is the part that makes it hold. A NEW OPERATOR DECLARES ITS KIND IN THE TYPE, NOT IN A COMMENT.**
The coefficient carries the distinction -- `Coefficient1p{OrdinaryKind}` against `Coefficient1p{ReducedKind}` -- and the
contraction is defined only for a matching pair, so pairing the wrong things raises a `MethodError` naming both types
instead of returning a wrong number. Do not add a coefficient type that omits the kind, and do not "fix" a mismatch by
inserting a compensating factor at the call site: that is the drift this rule exists to stop.

**Why documentation alone was tried and failed.** A careful note ALREADY existed at `module-Hfs.jl:370` describing this
very factor -- and the defect shipped in that same module, putting an uncorrected `sqrt(2j_a+1)` into every M1, E2 and
M3 hyperfine amplitude. Four further modules were bitten independently. The reason is worth stating plainly, because it
is the real lesson: **two quantities that differ by a factor obvious to whoever derived them get written with the same
symbol, and the factor then migrates to wherever it is convenient.** The similarity of the notation, not the difficulty
of the physics, is what has made every re-implementation of this machinery hard. A type cannot be talked into
forgetting; a convention in prose can.

**Fourth, if a compensation must live inside an integral, say so AT the integral.** `module-LandeZeeman.jl:227` uses
`coeff.T` bare with a rank-1 operator and is CORRECT, because `InteractionStrength.zeeman_n1` carries the
`1/sqrt(2j_a+1)` within itself. Nothing at the call site shows this. Any such integral must carry the factor in its own
docstring, or the next migration will edit the visible sites and silently break the invisible one.

**The worked example, kept current rather than described here:** the migration inventory in the module docstring of
`src/module-SpinAngularNew.jl` lists all eighteen sites by file and line, and separates those whose factor is visible
from those whose is not. Re-derive it rather than trusting the list -- `grep -n "coeff.T" src/*.jl` is the whole method,
and a coefficient used with no visible factor is a question, not an answer.

## Migrating the spin-angular convention: one commit, or none (Rule 19)

Written 25-Aug-2026, with the migration itself deliberately POSTPONED. The plan is recorded here rather than in a
message because the danger is not that it is hard but that it is done PARTIALLY, and a half-done convention change
produces wrong numbers with nothing failing.

**Why it cannot be staged.** Eighteen sites across nine modules compensate for the present convention, in two OPPOSITE
directions (see the inventory in the module docstring of `src/module-SpinAngularNew.jl`). Any intermediate state has
some callers on the old convention and some on the new, and every result computed in between is wrong. This is the one
place where the one-module-per-task rule is suspended DELIBERATELY, and it must be said out loud in the commit rather
than worked around.

**The order, and none of it is optional:**

1. **Capture a physical baseline FIRST** -- the full suite plus one real case per affected module: a hyperfine
   constant, an isotope shift, a Lande factor, a crystal-field splitting, with the numbers written down. Without this,
   "nothing changed" cannot be checked, and the suite alone will not do it: the 54 tests do not cover every affected
   module.
2. **Show the complete diff before applying it.** Eighteen near-identical deletions of a `sqrt(...)` factor is exactly
   the shape an eye slides over. The bulk-change rule applies with full force.
3. **Handle the INVISIBLE compensations explicitly, and expect them to be the failure.** `module-LandeZeeman.jl:227`
   uses `coeff.T` bare with a rank-1 operator and is CORRECT, because `InteractionStrength.zeeman_n1` carries the
   `1/sqrt(2j_a+1)` inside the integral. It needs NO call-site edit and DOES need the integral changed -- the exact
   inverse of every other site. A further twenty-one bare uses of `coeff.T` each need the same reading of what the
   integral beside them returns.
4. **Flip the module last**, only once the callers are consistent.
5. **Verify against the BASELINE, not against a green suite.** A suite pass here is necessary and nowhere near
   sufficient.

**Two decisions to take before starting, because they change the work:** whether the memoisation caches
(`SHELL_A_CACHE`, `SHELL_W_CACHE`, `PARTNER_CACHE`) are scoped or made task-local at the same time -- they are
unbounded, session-lived and not thread-safe -- and whether `SpinAngular` is retired immediately or kept beside the new
module for a release. Keeping both is safest and doubles the surface.

## The four working lists (Rule 20)

Four lists carry everything that is outstanding. **None of them lives in this file, and none ever will**: they
name collaborators, unpublished work and correspondence, which Rule 14 keeps out of a published document. They
live in the assistant's own memory directory. What belongs here is the MECHANISM — what each list is for, how it
is numbered, and how something leaves it — so that the rules survive even when the contents cannot be shown.

| list | holds | ordered by | numbers |
|---|---|---|---|
| **priority** | defects and open work in `src/` | urgency, within sections | plain integers |
| **application** | physics driven from OUTSIDE — someone is waiting, a paper is in preparation, a reply is owed | the strength of the outside obligation | `A1, A2, …` |
| **competition** | capability other codes have and JAC does not | promise against effort | `C1, C2, …` |
| **challenges** | held open: we do not want to discard it, but WHEN, WHETHER and HOW it returns is unknown | nothing — it is a holding place | its ORIGINAL priority number |

**The prefixes exist so that a number always means one thing.** "Item 8" is the priority list, "A3" the
application list, "C11" the competition list. A challenge keeps the number it had as a priority item, so that
something which comes back comes back as itself and nobody has to work out whether it is the same thing twice.

**NUMBERS ARE PERMANENT IDENTIFIERS.** Closing an item deletes its text and RETIRES its number, leaving a gap.
The gap is correct and must not be tidied away: renumbering makes "item 8" mean something different from one day
to the next, so that neither side can refer to an item reliably. This holds across all four lists.

### What an item says

A revisited item states four things, because a line naming a defect without the other three cannot be ranked
against anything:

    PROBLEM   what is wrong, concretely enough to recognise it in the code
    EFFORT    XS = an hour or two;  S = days;  M = one to three weeks;  L = one to three months;  XL = longer
    GAIN      what changes if it is fixed — and "no number changes" is a legitimate and important answer
    VERIFIED  when its central claim was last checked against `src/`, and what was checked

An item with no EFFORT line has never been sized; one with no VERIFIED line has never been checked. Both are
worth seeing at a glance, so neither is filled in speculatively.

### The five ways an item can leave

1. **CLOSED** — the work was done. The commit is the record.
2. **CLOSED BY A FACT** — a premise proved false, or the work turned out to be already done. These cannot return,
   and the fact that ended them is written down so the question is not reopened.
3. **NOT WORTH DOING** — a real problem whose gain does not justify its effort. This is a legitimate end state,
   equal to the others, and the reasoning is recorded so that the next reader does not re-derive it.
4. **MOVED TO CHALLENGES** — kept, with the condition for its return written down.
5. **RETIRED AT THE MAINTAINER'S REQUEST** — his decision, and it is not revisited.

### The rule that all of this exists for: NOTHING IS PUT BACK "FOR TRIAGE"

*Triage* — from French *trier*, to sort — is the sorting of cases by urgency so that limited effort goes where it
helps most. It is deliberately a judgement about **which** cases deserve attention, made without treating any of
them. Applied to a list it means a pass that only asks: is this real, is it urgent, can it go?

The failure it names is this: **something already decided against is put back on the list so that it can be judged
again**. That quietly reverses the decision, because anything on the priority list is by default something to be
done. It costs nothing to write and costs the next reader real time.

It has happened here. A module was deleted from the priority list on 18-Aug-2026 with the reason "a decision NOT
to work on something is not an open item", and was carried again a short while later "so the decision is visible";
another item recorded in its own text that it had been "restored for triage". **The challenges list is the remedy**:
a decision that is worth remembering is visible THERE, with its condition, and something returns to the priority
list when that condition is MET — not when a reader wonders about it again.

### And the failure that is not this one

Between 27 and 28-Aug-2026, five items sent work at things already done or explicitly forbidden. **None was a
re-proposal; all five were entries that had gone stale while the code moved underneath them.** That is what the
VERIFIED line is for. An item whose text cannot be trusted is worse than no item, because it is acted upon.


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

### /testExamples

**Sample 60 example branches at random, run each one alone, and JUDGE the result — then put what fails on the
priority list.** Meant to be run at night. It takes roughly 45 minutes of compute and needs no supervision.

It exists because the test suite and the examples fail differently. The suite tests what somebody thought to
test; the examples are where the physics is actually used, and a branch nobody has run in six months drifts
against the code silently. A first sweep of 200 branches on 28/29-Aug-2026 found two live regressions, nine
wrong grids and two one-character typos, none of which any test could see.

**1. Choose the branches.** A branch is one `if true` / `if false` / `elseif true` / `elseif false` guard in
`examples/*.jl`; there are about 500. Draw 60 AT RANDOM, but from the LEAST-RECENTLY-VISITED first: keep the
list of what each run sampled, and prefer branches not drawn before. Pure random sampling revisits some and
never reaches others; drawing from the unvisited pool reaches every branch in about nine runs. Record the seed.

**2. Run each one ALONE.** Copy the whole example file, set that one guard to `true` and EVERY other guard in
the file to `false`, so exactly one branch executes. Four harness details, each learned by getting it wrong:

- **Prepend the packages.** Only about a quarter of the example files carry a `using` line; the rest assume a
  live session. Prepend `using JenaAtomicCalculator`, `const JAC = JenaAtomicCalculator`, and
  `using Printf, LinearAlgebra, Distributed`, with `JLD2` and `SymEngine` inside a `try`. Without `Printf`
  a batch of examples fails on `@printf` and it looks like a defect.
- **Julia resolves `include` relative to the SCRIPT's directory, not the working directory.** One example
  includes a sibling file; copy that file next to the generated scripts, and run with the working directory
  set there so relative reads and `zzz-*.sum` output land in scratch rather than in the repository.
- **Use `timeout -k`.** Plain `timeout` sends only SIGTERM, which a branch inside a long numerical call can
  ignore; one escaped for eight hours.
- **Two at a time**, no more — each Julia process holds a couple of gigabytes.

**3. VERIFY THE HARNESS BEFORE BELIEVING IT.** Run five branches first and read every failure. If a failure is
the harness rather than the code, fix it and start again. The first attempt at the 28-Aug sweep reported 21
failures of which more than half were its own; a list half made of noise is worse than no list.

**4. Judge each result — this is the part that makes the command worth running.** Completing is NOT passing.
For every branch, read what the branch itself claims and compare:

- **Its own recorded numbers.** Most branches carry a REPORT block or a `Last successful` note quoting values.
  Do they still come out? A branch that completes while its own quoted number has moved by 30 % is a FINDING,
  and it is invisible to any "does it run" check. **Two traps here, both of which produced a false finding on the
  first run: (a) MOST BRANCHES SEND THEIR RESULTS TO A SUMMARY FILE, not to stdout — a number absent from the log
  is usually in the `.sum` the branch opened, so search both; (b) A LAST-DIGIT DIFFERENCE IS NOT DRIFT —
  5.28415 against 5.28416 is the same number rebuilt, and only a change in the leading figures is a finding.**
  Expect this check to reach only a minority of branches: on the first run just 4 of 38 completed branches quoted
  numbers in a mechanically comparable form. **Say how many were actually checked; a branch that merely ran is
  NOT verified.**
- **Physical sense.** Right order of magnitude, correct quantum numbers, partial sums matching totals, gauge
  agreement in the expected range, rates that are not zero and not 1e30, no negative populations.
- **Its own warnings.** A run that completes while printing NOT CONVERGED, or a grid warning, has not passed.
- **What the date claims.** A branch marked `Last successful: <date>` that now fails or has drifted is a
  REGRESSION and ranks far above a branch marked `unknown`, which has simply never worked. Say which.

**5. Report in six classes, and do not blur them.**
`VERIFIED` (ran and its numbers stand) · `DRIFTED` (ran, numbers moved — quote both) · `FAILED` (raised) ·
`TIMED OUT` (not a defect, say so) · `NOT A DEFECT` (a branch needing an earlier branch of the same file, or a
guard correctly refusing a wrong grid — the code was right and the EXAMPLE is wrong; still worth an item, but
say plainly that `src/` behaved) · **`RAISED AS DESIGNED`**.

**That last class is the one a does-it-run sweep gets backwards, so read the branch's own text before judging it.**
Some branches are GUARD TESTS: they ask for something impossible on purpose, to check that the code refuses
instead of returning a plausible number. For those, raising is the PASS and completing would be the defect. The
28-Aug sweep marked `example-Pc.jl` branch b a failure when its own comment says "a guard test, not a physics
branch" and sets `width = 0` deliberately. Judge such a branch by whether the error it raises is the error it
says it expects — and if it dies on a DIFFERENT error first, that is a real finding, because the earlier fault
is masking the guard the branch exists to demonstrate. That is exactly what had happened there.

**6. Add DRIFTED, FAILED and example-is-wrong results to the priority list**, one item each, with the file and
line, the error or the two numbers, and **one sentence saying what made you doubt it** — that sentence is the
point of the exercise and is what a later reader will act on. Do NOT add timeouts, and do NOT add branches that
merely need an earlier branch to have run.

**7. Never re-date a branch from this command.** Rule 7 asks for verified physical consistency, and a sampled
overnight run is not that. `/testExamples` reports and files items; dating stays with `/test Mx` and with a
human deciding.

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
3. Fix comments, per Rule 15.6-15.8 AS AMENDED 28-Aug-2026: a standalone `##` becomes `#`, but an APPENDED
   `##` on a code line STAYS -- and is instead read, and DELETED outright if it is not instructive; pure
   hash-only separators are KEPT, and so is a separator carrying a SHORT section label such as
   `####### Radiative Recombination ####`; only a separator carrying a sentence or enclosing a comment block
   goes; empty `#` lines go. Then **re-flow docstring prose and
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
One of the four lists described in Rule 20; the others are reached by `/application`, `/challenges` and by asking
for the competition list by name.
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

### /challenges
**Display the challenges list — what is held open, and on what condition it comes back.**

The fourth list, beside `/priority` (defects in `src/`), the application list (physics driven from outside) and
the competition list (capability other codes have). Like them it lives in the assistant's memory directory, not
here, and this command describes only the mechanism.

**What it is good for.** Some work is neither to be done nor to be forgotten: a measured speed-up that falls below
the standing bar, a model whose formula could not be confirmed, an import nobody yet has files for, a module
postponed so long that nobody remembers what it would cost. Left on the priority list these read as tasks, and
someone eventually starts them. Deleted outright, the reasoning is lost and the same idea is proposed again in six
months. The challenges list is where such a thing is kept ON PURPOSE, with three questions answered — **WHEN** it
would come back, **WHETHER** it should, and **HOW** much it would be — so that the next reader inherits the
judgement instead of repeating it.

1. Read the challenges list and print it **in full**, with each item's ORIGINAL priority number.
2. For each, say whether its stated CONDITION has been met. That is the only question this command asks.
3. **Move an item back to the priority list only when its condition is met**, and say which condition and what
   changed. It returns under its own old number.
4. **Never move an item back "for triage"** — see Rule 20. Re-examining something already settled, with no new
   evidence, is the failure this list exists to prevent.
5. If a condition has become impossible, the item is closed and its number retired, with the reason recorded.

The command only reads, prints and checks conditions. It never starts work on an item.


### /summarize
**End-of-session structured summary.**
Lists: modules changed (and what kind of change), example branches tested (with dates),
open items (blank dates, known bugs, planned next steps).
Important findings are saved to the memory system.

### /preparetravel
**Prepare both machines before leaving. Steps are labelled with the machine they run on.**

The aim is that the laptop carries everything needed and the desktop is left in a state that can be returned to
without archaeology. The whole reason this command exists: work left uncommitted on a machine you are flying away
from cost an evening once, and the fix is thirty seconds.

**ON THE DESKTOP** — leave nothing behind:
1. `git status`. If it is not clean, decide file by file. Anything worth keeping is committed and pushed; anything
   unfinished goes to a branch rather than staying loose:
   `git checkout -b wip/desktop-<date> && git add -A && git commit -m "WIP" && git push -u origin HEAD`
2. `git push` whatever is on the current branch, so the laptop can reach it.
3. Note what is still RUNNING. Long jobs continue while you travel only if they were started under `tmux` or
   `screen`; a job started in a plain ssh window dies when the link drops.

**ON THE LAPTOP** — take copies:
4. `git pull` on every branch you may want to touch while away.
5. `jac-fetch` — brings `apps/` across. `jac-results` — brings `work/` output if you want to read it en route.
6. `jac-memory-pull` — brings the assistant's memory across, so nothing has to be re-explained.
7. **Test the VPN before leaving.** The desktop sits at a private address: it is reachable from the institute
   network, and NOT from a hotel without the VPN. That is much better discovered at home than in a departure lounge.

**What NOT to copy:** `Manifest.toml` is machine- and Julia-specific. If an environment misbehaves on either
machine, the fix is `rm Manifest.toml` and `Pkg.instantiate()`, never a copy from the other side.


### /returntravel
**Come home: fold the trip's work back and start again on the desktop.**

**ON THE LAPTOP** — hand everything over:
1. `git status`, then commit and push everything, including any branch made while away.
2. `jac-send` — returns `apps/` with whatever changed while travelling.
3. `jac-memory-push` — returns the memory, so the desktop sessions know what happened on the trip.

**ON THE DESKTOP** — take it up again:
4. `git pull`, and merge any travel branch into the line it belongs to.
5. `rm Manifest.toml && julia --project=. -e 'using Pkg; Pkg.instantiate()'` if Julia has moved on either side.
6. Run the suite from `test/` (Rule 11). It is the proof the desktop is ready, not a formality.

**Which machine wins.** For `apps/` and for the memory the rule is simply: **the machine you were working on wins**,
and these two commands push in that direction. `jac-fetch`/`jac-send` and the memory pair copy only what is newer or
missing and delete nothing, so a mistaken direction is recoverable -- but running them the wrong way round can carry
a stale file over a fresh one, so follow the labels.

**Starting a SECOND session.** Both sessions work in the **same directory, on the same branch**. No worktrees, no
second branch, no merging between sessions, and nothing for the maintainer to remember or run.

What makes that safe is a discipline that belongs to the SESSIONS, not to him. Each session:

- **never uses `git add -A`, and never commits without an explicit pathspec** -- `git commit -F <msg> -- <only the
  files this session changed>`. A commit without a pathspec on a shared tree sweeps up the other session's staged
  work; this has happened, twice in one week, and the pathspec is what stopped it;
- **never runs `git checkout`, `git stash`, `git pull` or `git merge` without saying so first**, since all four
  change files under the other session;
- **states at the start which modules are its own**, and stays off the others'. A module assigned to neither session
  belongs to neither;
- **does not run the test suite while the other session is running it** -- they write the same files under `test/`.

If a session cannot say which files are its own, it must not commit.

**PARTITION BY DEPENDENCY, NOT BY TODAY'S FILE LIST.** Naming the modules each session owns is necessary and not
sufficient: a session working a THEME will reach every module that theme touches, days before it opens them. On
27-Aug-2026 one session migrated the spin-angular convention through 25 modules in an afternoon, while the other
went looking for "independent" work and picked `ImpactExcitation` -- untouched at that moment, and three
`SpinAngular` references away from being next. The usable test is one grep: a module with ZERO references to the
other session's theme is safe, one with any is not, however quiet it looks. State the theme, not just the files.

**BEFORE `git commit`, THREE CHECKS, because the tree moves under you.**

- `git status --porcelain` and stage by explicit pathspec, never `-A`, `.` or `-u`;
- `git diff --cached <each file>` must contain ONLY your own hunks. This is what proves you are not reverting the
  other session's work in a file you both edited. On 27-Aug `module-Hamiltonian.jl` changed hands four times,
  including a revert and a re-land, while an edit of the other session's sat in the working copy the whole time;
  the staged diff showing two hunks and nothing else is what made that safe to commit;
- **re-read the code your commit MESSAGE describes.** A message drafted ten minutes ago can be false by the time
  it is written. On 27-Aug a commit explaining a `sqrt(2j+1)` at a call site was drafted, and the other session
  removed that very factor before it was committed; the message had to be amended before pushing. A message is a
  claim about the code as it stands at HEAD, not as it stood when the work was done.

**A TEST RUN MEANS NOTHING UNLESS THE TREE STOOD STILL FOR ALL OF IT.** Not running the suite at the same time is
only half the rule: an EDIT by the other session during the run is just as corrupting, and silent. Julia compiles
the package once at startup, so a run that began before an edit may or may not have seen it, and afterwards there
is no way to tell. Note the modification times of every file the other session has dirty before starting, check
them again at the end, and report the window: "suite 54/54, their files unchanged at HH:MM:SS across a run from
HH:MM to HH:MM". If the times moved, the run is not evidence and must be repeated. Three suite runs were wasted
this way on 27-Aug before the rule was written down.

**PREFER A TARGETED CHECK TO THE FULL SUITE while the other session is active.** The suite needs a still tree for
three minutes; calling the changed functions directly and comparing numbers needs seconds and cannot be
contaminated. Where a change is claimed to move no number, that claim is better proved by a direct before/after
on the affected functions -- and where a rebuild alone perturbs the result (Rule 12's neighbour, the build
sensitivity), the ONLY sound proof is a CONTROL: put the original file back, rebuild, run, and show it reproduces
the new numbers. That control was needed twice on 27-Aug and settled both cases in minutes.

**THE MEMORY DIRECTORY AND THE WORKING LISTS ARE SHARED TOO, and they have no merge.** Two sessions editing the
priority list or `MEMORY.md` will silently overwrite each other: on 27-Aug an index line was clobbered and had to
be restored, and an item appeared in the list between one session reading it and writing it back. So: re-read the
file immediately before editing it, edit by ANCHORED replacement of a known string rather than by line number,
never rewrite a whole shared file from a copy held in memory, and check afterwards that the other session's
entries are still there. If two sessions must both record something, each adds its own item rather than
renumbering or reflowing the other's.

**WHEN THE THEMES COLLIDE, SAY SO AND STOP RATHER THAN WORK AROUND IT.** A session that cannot find work outside
the other's theme should say that plainly and wait, or ask for a different assignment. Half a day of hunting for
something "independent" in a codebase one migration is sweeping produces proposals for work already done, and
that is worse than idleness: on 27-Aug four such proposals in a row turned out to be already fixed or explicitly
forbidden. Waiting is a legitimate answer; inventing scope is not.


### recall commands
Print the command list above (this section). No steps performed.

## Reference file update protocol

When a code change legitimately alters numerical output (a physics fix or new feature):
1. The new `.sum` output must be verified against an independent source before updating.
2. Updating a `test/approved/*.sum` reference file is a deliberate editorial act — never do it
   automatically or silently.
3. A purely structural change (rename, refactor) must pass all tests unchanged. If a test fails
   after such a change, something went wrong.
