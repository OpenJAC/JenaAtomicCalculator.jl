#!/usr/bin/env python3
"""Run a RANDOM SAMPLE of example branches, each one ALONE, and collect what each printed.

    tools/sweepExamples.py --seed 20102026 --draw 150 --outdir <scratch>        # the full tour
    tools/sweepExamples.py --seed 20102026 --draw 150 --outdir <scratch> --only 5   # verify the harness first

WHY THIS FILE EXISTS.  `/testExamples` had been rebuilt from its description on every tour, and the four harness
details that description lists were each learned by getting them wrong once.  They are encoded here instead:

  1. **THE PACKAGES ARE PREPENDED.**  Only about a quarter of the example files carry a `using` line;  the rest
     assume a live session.  Without `Printf` a batch fails on `@printf` and it reads as a defect in the code.
  2. **JULIA RESOLVES `include` RELATIVE TO THE SCRIPT'S DIRECTORY, not the working directory**, so a file that
     includes a sibling needs that sibling copied next to the generated script.
  3. **EVERY BRANCH GETS ITS OWN DIRECTORY.**  A branch's `.sum` is opened by name, so two concurrent runs of
     different branches would otherwise write the same file and the drift check would read the wrong one.
  4. **`timeout -k`, NOT PLAIN `timeout`.**  SIGTERM alone is ignored by a branch inside a long numerical call;
     one escaped for eight hours on the 28-Aug-2026 tour.

WHAT IT DOES NOT DO:  it does not judge.  It records exit code, wall time and the log, and leaves the verdict to
a reader -- `tools/checkExampleDrift.py` for whether a branch's own quoted numbers came back, and a human for
whether the physics is sensible.  Completing is NOT passing, and no column here should be read as if it were.
"""

import argparse, csv, os, random, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

GUARD  = re.compile(r"^(if|elseif)\s+(true|false)\s*$")
INCL   = re.compile(r"""include\(\s*["']([^"']+)["']\s*\)""")
HEADER = """using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra, Distributed
try  using JLD2       catch;  end
try  using SymEngine   catch;  end

"""


def guards(exampleDir):
    """Every branch guard in the tree, as (path, lineNumber) pairs in file then line order."""
    out = []
    for name in sorted(os.listdir(exampleDir)):
        if not (name.startswith("example-") and name.endswith(".jl")):
            continue
        path = os.path.join(exampleDir, name)
        with open(path, encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh, start=1):
                if GUARD.match(line.rstrip("\n")):
                    out.append((path, i))
    return out


def writeScript(path, line, runDir):
    """Copy the example with THIS guard true and every other guard in the file false."""
    os.makedirs(runDir, exist_ok=True)
    with open(path, encoding="utf-8", errors="replace") as fh:
        src = fh.readlines()
    for i in range(len(src)):
        m = GUARD.match(src[i].rstrip("\n"))
        if not m:
            continue
        want = "true" if i + 1 == line else "false"
        src[i] = re.sub(r"(true|false)\s*$", want, src[i].rstrip("\n")) + "\n"
    script = os.path.join(runDir, os.path.basename(path))
    with open(script, "w", encoding="utf-8") as fh:
        fh.write(HEADER)
        fh.writelines(src)
    # rule 2 above:  a sibling that the example includes has to sit beside the generated script
    for rel in INCL.findall("".join(src)):
        sib = os.path.join(os.path.dirname(path), os.path.basename(rel))
        if os.path.isfile(sib):
            shutil.copy(sib, os.path.join(runDir, os.path.basename(rel)))
    return script


def runOne(job):
    path, line, outdir, jacRoot, limit = job
    tag    = f"{os.path.basename(path)[:-3]}-L{line}"
    runDir = os.path.join(outdir, tag)
    script = writeScript(path, line, runDir)
    log    = os.path.join(runDir, "run.log")
    t0     = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        rc = subprocess.call(["timeout", "-k", "20", str(limit), "julia",
                              f"--project={jacRoot}", os.path.basename(script)],
                             cwd=runDir, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    with open(log, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    # 124 is timeout's own code;  137 is the SIGKILL that -k delivers afterwards
    klass = "TIMED OUT" if rc in (124, 137) else ("ran" if rc == 0 else "NONZERO EXIT")
    err   = ""
    for pat in ("ERROR: LoadError: ", "ERROR: "):
        k = text.find(pat)
        if k >= 0:
            err = text[k + len(pat):].split("\n")[0][:200]
            break
    warn = sum(text.count(w) for w in ("NOT CONVERGED", "not converged", "BOX TOO SMALL",
                                       "NOT SATURATED", "enlarge box"))
    print(f"   {klass:<12} {tag:<28} {dt:7.1f} s   rc={rc}" + (f"   {err[:70]}" if err else ""),
          flush=True)
    return dict(example=os.path.basename(path), line=line, tag=tag, klass=klass, rc=rc,
                seconds=round(dt, 1), warnings=warn, error=err, log=log)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed",    type=int, required=True)
    ap.add_argument("--draw",    type=int, default=150)
    ap.add_argument("--outdir",  required=True)
    ap.add_argument("--only",    type=int, default=0, help="run just the first N of the draw")
    ap.add_argument("--skip",    type=int, default=0, help="skip the first N of the draw")
    ap.add_argument("--timeout", type=int, default=500)
    ap.add_argument("--jobs",    type=int, default=2)
    ap.add_argument("--root",    default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    a  = ap.parse_args()

    allGuards = guards(os.path.join(a.root, "examples"))
    drawn     = random.Random(a.seed).sample(allGuards, min(a.draw, len(allGuards)))
    chosen    = drawn[a.skip:(a.skip + a.only) if a.only else len(drawn)]
    os.makedirs(a.outdir, exist_ok=True)

    print(f"seed {a.seed};  {len(allGuards)} guards in the tree;  {len(drawn)} drawn uniformly, "
          f"running {len(chosen)} of them (skip {a.skip}), {a.jobs} at a time, {a.timeout} s cap\n", flush=True)
    jobs = [(p, l, a.outdir, a.root, a.timeout) for p, l in chosen]
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        rows = list(pool.map(runOne, jobs))

    out = os.path.join(a.outdir, f"results-seed{a.seed}-skip{a.skip}.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader();   w.writerows(rows)
    n = {}
    for r in rows:
        n[r["klass"]] = n.get(r["klass"], 0) + 1
    print(f"\n{out}")
    print("   " + ",  ".join(f"{k}: {v}" for k, v in sorted(n.items())))
    print("   COMPLETING IS NOT PASSING -- the verdict still has to be read off each log.")


if __name__ == "__main__":
    main()
