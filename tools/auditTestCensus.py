#!/usr/bin/env python3
"""Measure which src/ modules the test suite reaches, and how.  READ-ONLY -- it changes nothing.

    usage:   python3 tools/auditTestCensus.py [--cov <dir>] [--all]

WHY A SEPARATE FILE FROM `auditHygiene.py`.  That one measures each source file on its own; this
question is per MODULE and has to read `test/runtests.jl` and `module-TestFrames*.jl` as well, so it
cannot be a per-file row there.  The register line of 07-Nov-2026 says the calibration for this
detector sits at the head of `auditHygiene.py`;  it does not, and this file is where it belongs.

THE QUESTION IS NOT ONE QUESTION BUT THREE, and conflating them is what makes every naive count
wrong.  A module may:

    (a) have a test carrying ITS OWN NAME        -- `testModule_Einstein` for `Einstein`
    (b) be ENTERED by the suite without one      -- `Basics` is entered by every test there is
    (c) be neither                               -- and only (c) is a gap

A count of (a) alone reads 35 of 108 modules as untested, where the modules nothing executes are
far fewer.  Pass `--cov` with a directory of `.cov` files from `julia --code-coverage=user` to
decide (b) by MEASUREMENT rather than by reference-counting;  without it, (b) is inferred, and the
inference is an over-estimate -- a reference may sit behind a guard that no test turns on.
`InteractionStrengthQED` had five call sites in `src/` and had never been executed by anything,
because every one of them sits behind `settings.qedModel != NoneQed()`.

Four calibrations, each paid for with a wrong answer:

    1. A TEST MAY BE CALLED DIRECTLY FROM `runtests.jl`, not through `TestFrames`.  `RacahAlgebra`
       is tested by four such calls and reads as untested if only `TestFrames.*` is collected.
       The priority list recorded this on 28-Aug-2026 -- "the old item 32 was stale from birth
       because the detector could not see a test called directly from `runtests.jl`" -- and it was
       re-derived, the same way, on 08-Nov-2026.  Hence: collect every `@test <name>()`.
    2. A MODULE MAY BE EXERCISED UNDER ANOTHER NAMESPACE.  `module-ManyElectronAZ.jl` defines
       `ManyElectron.matrixElement_*`, which the suite tests under the OTHER module's name, so the
       file has no test of its own and no gap either.  Hence: a module that defines nothing in its
       own namespace is reported separately.
    3. A TEST MAY NAME A TYPE RATHER THAN THE MODULE THAT HOLDS IT.  Five of `Plasma`'s ten schemes
       are tested -- `AverageAtom`, three `SuperConfiguration*`, `KineticChargeBalance` -- and all
       five are `Plasma.AbstractPlasmaScheme`, so "Plasma has no test" is false in the way that
       matters.  Hence: a test name is matched against the module's own exported TYPE names too.
    4. `module-X-inc.jl` IS NOT A MODULE.  `module-LSjj-inc.jl` is an include file of `LSjj`, which
       IS tested;  counting it as a module put the gap one too high on 08-Nov-2026.
    5. A `-` IN A `.cov` FILE MEANS NO COUNTER, NOT A NON-EXECUTABLE LINE, and a file in which
       nothing ran has NO `.cov` FILE AT ALL -- which is the strongest signal there is, not missing
       data.  Julia compiles lazily, so an uncalled method's whole body reads exactly like a block
       of comments.  See `coveredLines`, which also has the converse trap.
    6. THE COVERED COLUMN IS TWO NUMBERS AND NOT A PERCENTAGE, because only what is COMPILED is
       instrumented:  ran/compiled is near 1 by construction, and `Math` read 100 % on one line of
       129 -- `polylogExp`, which the suite enters 144249 times and never checks the value of.
"""

import os, re, sys, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, "src")
TEST = os.path.join(ROOT, "test", "runtests.jl")

## THE REFLECTION TESTS, which enter EVERY module and test the physics of none.
## `testMethod_SettingsCopyConstructors` walks the submodules by reflection and calls every keyword
## copy-constructor it finds -- deliberately, because a hand-written list would go stale the moment a
## module is added -- and `testStructConstructors` does the same for the struct constructors and
## `Base.show`.  So EVERY module holding a `Settings` is entered by the suite, and reporting that as
## "the suite reaches it" is true and useless:  it put 20 of today's 34 gaps in the covered column
## while testing not one number in any of them.  They are counted apart, as SHALLOW.
REFLECTION = {"TestFrames.testMethod_SettingsCopyConstructors", "TestFrames.testStructConstructors"}

## The suffixes JAC appends to a physics name when it makes a type of it, so that a test naming the
## PHYSICS is matched to the type that carries it.  Kept short and explicit for the reason in
## calibration 3:  anything looser matches on a shared tail and hides a real gap.
SUFFIX = ["", "Scheme", "Settings", "Outcome", "Line", "Channel", "Computation", "Model"]


def stripped(path):
    """The file with docstrings and comments removed, so a MENTION is never read as a CALL.
    A module named only in a comment or a docstring is not reached by anything."""
    out, inD = [], False
    for l in open(path, encoding="utf-8", errors="replace"):
        s = l.lstrip()
        if not inD and s.startswith('"""'):
            inD = not ('"""' in s[3:]);  continue
        if inD:
            inD = '"""' not in l;  continue
        out.append(re.sub(r'#.*$', '', l))
    return "\n".join(out)


def modulesOf():
    """-> {module name: [its files]}.  An `-inc-` file and a bare `-inc` file both belong to the
    module they are named after;  neither is a module (calibration 4)."""
    mods = collections.defaultdict(list)
    for f in sorted(os.listdir(SRC)):
        if not (f.startswith("module-") and f.endswith(".jl")):  continue
        mods[re.sub(r'-inc(-.*)?$', '', f[len("module-"):-len(".jl")])].append(f)
    return mods


def testNames():
    """-> every name the suite asserts, as a set of underscore-separated tokens per test.
    Collected from `@test <anything>()`, NOT from `@test TestFrames.*` (calibration 1)."""
    txt = open(TEST, encoding="utf-8").read()
    txt = re.sub(r'#.*$', '', txt, flags=re.M)
    out = {}
    for m in re.finditer(r'@test\s+([A-Za-z_][\w.]*)\s*\(\s*\)', txt):
        full = m.group(1)
        name = full.split(".")[-1]
        body = name.split("_", 1)[1] if "_" in name else name
        toks = set(body.split("_"))
        ## THE QUALIFIER NAMES A MODULE TOO.  `@test RacahAlgebra.testSumRules()` says which module
        ## it tests in the part before the dot and nowhere else -- its own name carries no
        ## underscore and so no token.  Reading tokens alone reported RacahAlgebra as untested,
        ## which is calibration 1 above, made twice.
        if "." in full and full.split(".")[0] != "TestFrames":  toks.add(full.split(".")[0])
        out[full] = toks
    return out


def typesOf(files):
    """-> the type names a module defines, so a test naming a TYPE counts (calibration 3)."""
    names = set()
    for f in files:
        for m in re.finditer(r'^\s*(?:struct|abstract\s+type|mutable\s+struct)\s+([A-Za-z_]\w*)',
                             stripped(os.path.join(SRC, f)), re.M):
            names.add(m.group(1))
    return names


def ownNamespace(mod, files):
    """False when every function the module defines belongs to ANOTHER module (calibration 2)."""
    mine = other = 0
    for f in files:
        for m in re.finditer(r'^function\s+(?:([A-Za-z_]\w*)\.)?[A-Za-z_]', stripped(os.path.join(SRC, f)), re.M):
            if m.group(1) and m.group(1) != mod:  other += 1
            else:                                 mine  += 1
    return mine > 0 or other == 0


def coveredLines(covDir, files):
    """-> (executed, instrumented) line counts over a module's files, or None when no --cov was given.

    TWO THINGS ABOUT `.cov` FILES THAT INVERT THE READING, and both were got wrong on 08-Nov-2026.
    A `.cov` line is `       12 code` for twelve executions and `       0 code` for an executable line
    never reached -- but `       - code` means NO COUNTER EXISTS, and Julia compiles lazily, so a
    method that is never called is never compiled and its whole body comes back as `-`, exactly like a
    comment.  Counting only the `0` lines therefore reported `Math`, `Distribution`, `Beam` and
    `AtomicFeatures` as having no unentered function at all, when `Math` has ONE instrumented line in
    129.  And a file in which nothing ran at all gets NO `.cov` FILE, which read as missing data rather
    than as the strongest possible signal:  `module-ForPedestrians.jl`, 1539 lines and 20 functions,
    has none.  So a missing file counts as zero instrumented lines, not as unknown.

    The converse trap is `module-Basics.jl`:  134 lines holding 103 bodiless `function name end`
    declarations and three includes, so it has no executable line of its own and correctly has no
    `.cov`.  Its code lives in the `-inc-` files, which is why the counts are rolled up per MODULE
    over all of its files rather than read per file."""
    if not covDir:  return None
    hit = run = 0
    for f in files:
        for c in sorted(os.listdir(covDir)):
            if not (c.startswith(f + ".") and c.endswith(".cov")):  continue
            for l in open(os.path.join(covDir, c), encoding="utf-8", errors="replace"):
                tag = l[:9].strip()
                if tag.isdigit():
                    run += 1
                    if int(tag) > 0:  hit += 1
    return (hit, run)


def main():
    covDir = sys.argv[sys.argv.index("--cov")+1] if "--cov" in sys.argv else None
    mods   = modulesOf()
    tests  = testNames()
    rows   = []
    for mod, files in sorted(mods.items()):
        tps  = typesOf(files) | {mod}
        ## A TOKEN MATCHES A TYPE EXACTLY, OR WITH ONE OF JAC'S OWN TYPE SUFFIXES APPENDED:
        ## `testModule_AverageAtom` names `Plasma.AverageAtomScheme`, and five of Plasma's ten schemes
        ## are tested that way, so requiring bare equality reported Plasma as untouched.  But a plain
        ## SUBSTRING test is wrong in the opposite direction and worse, because it hides a gap:  it
        ## credited `testModule_PhotoExcitation` with reaching `BeamPhotoExcitation`, and
        ## `testModule_AutoIonization` with reaching `DoubleAutoIonization`, on nothing but a shared
        ## tail.  A token must name the type, not sit inside a longer one.
        hits = [t for t, toks in tests.items()
                if (toks & tps) or any(k + sfx in tps for k in toks for sfx in SUFFIX)]
        own  = [t for t in hits if t not in REFLECTION]
        ## The two reflection tests reach a module BY CONSTRUCTION, not by naming it:  one calls
        ## every `Settings` keyword copy-constructor it can find, the other every struct constructor
        ## and `Base.show`.  So the test is whether the module HAS such a type, and matching their
        ## names against it decides nothing -- which is why this is asserted here rather than left
        ## to `hits`.
        refl = sorted(REFLECTION) if tps - {mod} else []
        byName = [t for t in own if mod in tests[t]]
        cov  = coveredLines(covDir, files)
        rows.append((mod, files, byName, own, ownNamespace(mod, files), cov, refl))

    print("test census over %d modules in src/,  %d tests asserted by test/runtests.jl" % (len(mods), len(tests)))
    print("-" * 104)
    gapNamed = [r for r in rows if not r[2]]
    gapAny   = [r for r in rows if not r[3]]
    shallow  = [r for r in rows if not r[3] and r[6]]
    print("   modules with a test carrying their own name        %4d" % (len(mods) - len(gapNamed)))
    print("   modules with no test of their own name             %4d" % len(gapNamed))
    print("      of those, a test names one of their TYPES       %4d" % (len(gapNamed) - len(gapAny)))
    print("      of those, nothing in the suite names them       %4d" % len(gapAny))
    print("         of those, entered ONLY by the reflection tests %3d  (a constructor, no physics)"
          % len(shallow))
    print("   modules defining only another module's namespace   %4d" % sum(1 for r in rows if not r[4]))
    if covDir:
        never = [r for r in rows if r[5] and r[5][0] == 0]
        print("\n   MEASURED by coverage in %s" % covDir)
        print("   modules with NO line ever executed                 %4d" % len(never))
        for mod, files, byName, own, ownNs, cov, refl in sorted(never):
            print("      %-32s  %s" % (mod, "%d instrumented lines, none reached" % cov[1]
                                                  if cov[1] else "not one line instrumented: never compiled"))
    print("\n" + "-" * 104)
    hdr = "   %-32s %-9s %-11s %s" % ("module", "own name", "ran/compiled", "the test that reaches it")
    print(hdr);  print("-" * 104)
    for mod, files, byName, own, ownNs, cov, refl in rows:
        if byName and "--all" not in sys.argv:  continue
        ## THE TWO NUMBERS, NOT A PERCENTAGE.  Julia instruments only what it COMPILES, and it
        ## compiles only what is called, so ran/compiled is close to 1 by construction and a
        ## percentage flatters exactly where it matters:  `Math` read 100 % on ONE line of 129,
        ## the body of `polylogExp`, which the suite enters 144249 times and never checks.  The
        ## compiled count is the informative half and must stay visible.
        pc = "-" if not cov else ("%d/%d" % cov if cov[1] else "none")
        why = ", ".join(sorted(own)[:2]) if own else \
              ("constructor only (reflection)" if refl else
               "another namespace" if not ownNs else "")
        print("   %-32s %-9s %-11s %s" % (mod, "yes" if byName else "NO", pc, why))


main()
