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

## Reference file update protocol

When a code change legitimately alters numerical output (a physics fix or new feature):
1. The new `.sum` output must be verified against an independent source before updating.
2. Updating a `test/approved/*.sum` reference file is a deliberate editorial act — never do it
   automatically or silently.
3. A purely structural change (rename, refactor) must pass all tests unchanged. If a test fails
   after such a change, something went wrong.
