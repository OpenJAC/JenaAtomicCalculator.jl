#!/usr/bin/env python3
"""Check whether the numbers an example branch QUOTES IN ITS OWN COMMENTS still come out when it is run.

    tools/checkExampleDrift.py  examples/example-Ed.jl  8  /path/to/run.log  [extra.sum ...]

WHY THIS EXISTS -- priority item 78.  The `/testExamples` sweep's distinctive question is not "does the branch
run" (the test suite largely covers that) but "does a branch's own RECORDED number still come out".  That is the
one failure mode nothing else can see: a branch that still completes while its number has moved.

THE CHECK HAD NEVER WORKED, AND THE REASON IS WORTH STATING, because it is the whole design of this tool.  On
17-Oct-2026 a hand-rolled detector flagged 53 of 76 branches -- including `example-Ed.jl`, whose eight numbers
had been verified by hand the day before.  Calibrated afterwards on that same branch's output, by drawing random
values of the right magnitude and asking how often they "match" something the run printed:

      match within 1e-5   0.1 % of random values        <-- a real verdict
      match within 1e-3   3.4 - 4.1 %                   <-- weak;  about 1 in 25
      match within 1e-1    77 %                         <-- MEANINGLESS

53 of 76 is 70 %.  The detector's finding rate WAS its false-match rate.  With a few hundred numbers in a log,
a relative window of 10 % is nearly tiled, so almost any value finds a neighbour.

HENCE THE THREE RULES THIS TOOL IS BUILT ON:

  1. **THERE IS NO "DRIFTED" VERDICT, AND THERE CANNOT BE.**  A number that moved by 5 % and a number that was
     never a printed result at all -- a literature value, an input, a ratio, a date -- look IDENTICAL from here.
     So the only verdict is REPRODUCED, and everything else is NOT FOUND, which is a question and not a finding.
  2. **THE TOOL CALIBRATES ITSELF ON EVERY RUN** and prints the expected false count beside the real one.  A
     report that says "11 of 42 reproduced, expected by chance 0.0" means something;  one that says "19 drifted"
     does not.
  3. **A LAST-DIGIT DIFFERENCE IS NOT DRIFT.**  5.28415 against 5.28416 is the same number rebuilt, so the
     closest printed number and the relative difference are always shown and the reader judges the margin.
  4. **A RELATIVE WINDOW HIDES A MOVE IN A SMALL DIFFERENCE, so the ABSOLUTE difference is shown too.**  Found
     on the 18-Oct-2026 sweep, and it is the one way this tool can report REPRODUCED about a claim that has in
     fact moved.  `example-Ao.jl` branch at line 152 quotes an EOL total of -288.66499124 Ha and got
     -288.66401063 -- 3.4e-06 relative, so REPRODUCED -- but the claim's own subject is the EOL-AL GAIN, and
     that gain fell from -0.00393584 to -0.00295523 Ha, a quarter of itself.  Both numbers are near 288 and
     their DIFFERENCE is near 0.003, so no relative test on the totals can see it.  The absolute column shows
     9.8e-04 Ha and a reader asks the right question.  **Where a claim's subject is a DIFFERENCE of two large
     numbers, check the difference and not the totals** -- the same rule that says not to subtract large
     near-equal numbers, read from the other end.

FOUR MEASURED SOURCES OF NOISE, each handled below and each found by running the control:
  * **A PROSE EXPONENT.**  A table writes a column of mantissas and ONE exponent in the words -- "2.212102 /
    1.734163 / ... e5 barn".  Comparing a mantissa against a full value produced the original false findings, so
    a standalone `eN`, `e+NN`, `10^N` or `10^-N` is taken as a shared exponent for the bare numbers on its line.
  * **A PROSE RESTATEMENT.**  The same branch quotes `2.212102` in its table and `2.212` in the sentence below
    it.  A quoted number whose digits are a PREFIX of another quoted number's is the same value said twice, and
    only the longest form is kept.
  * **DATES AND COUNTS.**  `15-Oct-2026` yields `-2026`;  `3000 eV` is an input.  A recorded physics result
    carries a decimal point, so bare integers are not read at all.
  * **ANOTHER BRANCH'S OUTPUT.**  Only the `.sum` THIS example opens is read, found from its own
    `setDefaults("print summary: open", ...)` call -- never every `.sum` lying in a directory.
  * **A BRANCH'S HISTORY, WHICH IS THE WORST OF THEM.**  A comment block holds the current claim AND every
    superseded value, convergence table and literature target.  Checking the history is how a value that was
    deliberately replaced becomes a "finding", so where the `Last successful` line quotes numbers, only those
    are checked.  Two traps found by the controls: a block may carry a STALE `Last visit ... NOT successful`
    line above its live one, so `Last successful` wins and the latest wins;  and the claim is capped at a few
    lines, because on a block with no blank-comment separator the run-on swallowed a paper's whole target table.
  * **AN IDENTIFIER IS NOT A MEASUREMENT.**  `arXiv:2005.00414` parses as a number and was duly reported as a
    missing result.  So does a DOI.  Both are known only by what precedes them.

WHAT IT MEASURED ON 150 BRANCHES, 17-Oct-2026, so that the next sweep knows what to expect:

      67 of 150   quote no number comparable in this way at all -- the usual case, and NOT a finding
      83 of 150   carry a checkable claim;  32 of those reproduced with no question left over
     602 numbers  checked, 216 reproduced to better than 1e-5

AND THE ONE RESULT WORTH ACTING ON:  a branch whose `Last successful` line QUOTES ITS OWN NUMBERS leaves
**1.0** question per branch;  a branch whose numbers are buried in the prose history leaves **4.1**, because
there the whole block has to be read and its superseded values, literature targets and inputs come with it.
So the cheapest way to make a branch machine-checkable is one line of writing, not any more code here: put the
result ON the date line.  `examples/example-Dh.jl` branch (b) is the model --

    # Last successful:  15-Oct-2026   ## Coulomb 8.22865, Babushkin 8.22902 against the exact 8.22906 /s,

-- and it checks in four numbers with nothing left to read.

AND TEST IT ON A BRANCH KNOWN TO BE GOOD BEFORE BELIEVING A FINDING.  `examples/example-Ed.jl` (8 numbers
verified 16-Oct-2026) and `examples/example-Dh.jl` branch (b) make ready-made controls.  A detector that cannot
reproduce a known-good branch has no claim on an unknown one.
"""
import bisect, math, os, random, re, sys

# A recorded result carries a decimal point.  Bare integers are deliberately NOT matched: they are dates,
# counts and inputs, and they were a measured source of false findings.
NUM     = re.compile(r'[-+]?\d+\.\d+(?:[eE][-+]?\d+)?')
ANYNUM  = re.compile(r'[-+]?\d+\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?')
# a standalone exponent token:  e5   e+05   e-3   10^5   10^-3   x 10^12
# an identifier immediately before a number:  arXiv:2005.00414,  doi:10.1016/...,  (arXiv: ...)
IDENT   = re.compile(r'(?:arxiv|doi|isbn|arXiv)\s*:\s*(?:[-\w./]*)$', re.IGNORECASE)
LONEEXP = re.compile(r'(?:^|[\s(/,])(?:[eE]([-+]?\d{1,3})|(?:[x×]\s*)?10\^?\s*([-+]?\d{1,3}))(?=[\s).,;:]|$)')

REPRODUCED = 1.0e-5          # measured false-match rate 0.1 %
LASTDIGITS = 1.0e-3          # measured false-match rate ~4 %;  reported, never counted as a finding
CLAIMLINES = 4               # a `Last successful` line plus its continuations -- never a whole block


def quotedNumbers(path, guardLine, lookahead=80):
    """The numbers this branch quotes in its own comments, with a prose exponent applied and restatements dropped.

    A list of `(text, [candidate values])` is returned, longest-form first.
    """
    lines = open(path, errors='replace').read().split('\n')
    # `guardLine` is the 1-based line of the branch's own guard;  tolerate being handed that line or the one
    # after it, and in either case begin BELOW the guard rather than on it.
    start = guardLine
    if not (start < len(lines) and re.match(r'^\s*(if|elseif|else)\b', lines[start - 1])):
        if start >= 1 and re.match(r'^\s*(if|elseif|else)\b', lines[start]):
            start += 1
    block = []
    for l in lines[start:start + lookahead]:
        if re.match(r'^\s*(if|elseif|else)\b', l):
            break
        block.append(l)

    # A BRANCH'S COMMENT BLOCK HOLDS BOTH ITS CURRENT CLAIM AND ITS HISTORY, and checking the history is how a
    # superseded value becomes a "finding".  Rule 7 already marks the current claim: the `Last successful:` (or
    # `Last visit:`) line and its continuations are the dated statement of what this branch produces.  So where
    # that line carries numbers, THEY are what is checked and the remaining prose is context;  where it carries
    # none -- the common case -- the whole block is checked, as before.
    # A BLOCK MAY CARRY SEVERAL DATE LINES -- the control branch holds a stale `Last visit: ... NOT successful`
    # above its live `Last successful`, and taking the first found the stale one.  A `Last successful` therefore
    # beats a `Last visit`, and among equals the LAST one wins, history accumulating above the current claim.
    # AND A DATE LINE READING `unknown` IS NOT A CLAIM.  Measured on the 18-Oct-2026 sweep: two branches were
    # labelled dated claims they do not make -- `example-Dq.jl` says `Last successful:  unknown ...` and
    # `example-Cl.jl` says in its `Last visit` line that it CANNOT BE DATED, so in both cases the numbers
    # harvested came from the literature table below.  Neither produced a finding, but both wasted a reading.
    # AND ONLY `Last successful` IS A CLAIM.  Rule 7 as rewritten 17-Oct-2026 gives `Last visit` one meaning --
    # the branch RAN on that date, with no claim at all about the result -- so harvesting numbers from it
    # contradicts the rule;  `Last reproduced` is this tool's own output and quotes counts, not values.  Where
    # there is no `Last successful` the whole block is checked, as for any undated branch.
    claim, claimKind = [], None
    for pattern, kind in ((r'Last\s+successful\s*:', 'Last successful'),):
        hits = [i for i, l in enumerate(block)
                if re.search(pattern, l) and 'unknown' not in l.lower() and 'cannot be dated' not in l.lower()]
        if not hits:
            continue
        claimKind = kind
        # AND CAP IT.  Where a block carries no bare `#` separator the run-on would swallow the whole thing
        # and then be LABELLED a dated claim, which is worse than not finding one: on `example-Cb.jl` it ate a
        # paper's Dirac-Fock target table.  A dated claim is a line plus a continuation or two.
        for m in block[hits[-1]:hits[-1] + CLAIMLINES]:
            if '#' not in m or m[m.index('#'):].strip().strip('#').strip() == '':
                break
            claim.append(m)
        break

    # A PROSE LIST WRAPS, AND ITS EXPONENT SITS ON THE LAST LINE.  The control branch writes
    # "-- 2.212102 / 1.734163 /" on one line and "/ 1.316345 e5 barn" on the next, so a per-line exponent
    # reaches six of its eight numbers and misses two.  Exponents are therefore collected per contiguous
    # COMMENT PARAGRAPH and offered as extra candidates -- and the calibration is scaled by how many candidates
    # each number got, so widening the search cannot quietly inflate the verdict.
    def harvest(src_lines):
        paras, cur = [], []
        for l in src_lines:
            if '#' in l:
                cur.append(l[l.index('#'):])
            elif cur:
                paras.append(cur);  cur = []
        if cur:
            paras.append(cur)

        raw = []
        for para in paras:
            paraExp = set()
            for txt in para:
                for m in LONEEXP.finditer(txt):
                    paraExp.add(int(m.group(1) if m.group(1) is not None else m.group(2)))
            for txt in para:
                for t in NUM.finditer(txt):
                    s = t.group(0)
                    # AN IDENTIFIER IS NOT A MEASUREMENT.  `arXiv:2005.00414` parses as 2005.00414 and was
                    # reported as a missing result on `example-Cb.jl`;  so do DOIs.  Both are recognised by
                    # what PRECEDES them, since nothing about the digits themselves gives it away.
                    if IDENT.search(txt[:t.start()]):
                        continue
                    digits = s.lower().split('e')[0].replace('-', '').replace('+', '').replace('.', '').strip('0')
                    if len(digits) < 4:                   # too few digits to carry information
                        continue
                    try:
                        v = float(s)
                    except ValueError:
                        continue
                    if v == 0.0:
                        continue
                    cand = [v]
                    if 'e' not in s.lower():
                        for e in sorted(paraExp):
                            if v * 10.0 ** e not in cand:
                                cand.append(v * 10.0 ** e)        # the prose exponent, applied
                    raw.append((s, cand, digits))

        # A PROSE RESTATEMENT: `2.212` under a table that already said `2.212102` is the same value twice.
        raw.sort(key=lambda r: -len(r[2]))
        out, seen = [], []
        for s, cand, digits in raw:
            if any(d.startswith(digits) for d in seen):
                continue
            seen.append(digits)
            out.append((s, cand))
        return out

    dated = harvest(claim) if claim else []
    if dated:
        return (dated, f'the dated claim on the `{claimKind}` line')
    return (harvest(block), 'the whole comment block (its date line quotes no number)')


def summaryFiles(path):
    """The .sum files THIS example opens -- read from its own setDefaults call, not guessed."""
    txt = open(path, errors='replace').read()
    return re.findall(r'print summary: open"\s*,\s*"([^"]+)"', txt)


def printedNumbers(paths):
    """Every non-zero number the run printed, as a sorted list of absolute values."""
    vals = set()
    for p in paths:
        if not p or not os.path.exists(p):
            continue
        for t in ANYNUM.finditer(open(p, errors='replace').read()):
            try:
                v = float(t.group(0))
            except ValueError:
                continue
            if v != 0.0:
                vals.add(abs(v))
    return sorted(vals)


def closest(x, printed):
    """The printed number nearest `x`, and the relative difference;  `(None, None)` if nothing was printed.

    `printed` is sorted, so the neighbours are found by bisection rather than by a scan:  the calibration calls
    this a few thousand times per branch and a linear scan makes a bulk sweep needlessly slow.
    """
    if not printed:
        return (None, None)
    a = abs(x)
    i = bisect.bisect_left(printed, a)
    best, bestRel = None, None
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(printed):
            rel = abs(printed[j] - a) / a
            if bestRel is None or rel < bestRel:
                best, bestRel = printed[j], rel
    return best, bestRel


def falseMatchRate(printed, nCand=1.0, draws=2000, seed=20261017):
    """How often a RANDOM value of the right magnitude matches something printed -- this run's own calibration.

    `nCand` is the mean number of candidate scalings each quoted number is allowed, and the draw mimics that:
    a trial counts as a match if ANY of its `nCand` values matches, so offering more candidates RAISES the
    measured false rate instead of quietly improving the verdict.  A pair `(rate within REPRODUCED, rate within
    LASTDIGITS)` is returned.  Without this calibration a count of matches cannot be read at all, which is the
    defect that made the first attempt at this tool useless.
    """
    if not printed:
        return (0.0, 0.0)
    rng, c5, c3 = random.Random(seed), 0, 0
    k = max(1, int(round(nCand)))
    for _ in range(draws):
        hit5 = hit3 = False
        for _ in range(k):
            base = rng.choice(printed)
            ex   = math.floor(math.log10(base)) if base > 0 else 0
            x    = float('%.5e' % (rng.uniform(1.0, 9.999999) * 10.0 ** ex))
            if x == 0.0:
                continue
            _, rel = closest(x, printed)
            if rel is None:
                continue
            hit5 |= rel < REPRODUCED
            hit3 |= rel < LASTDIGITS
        c5 += hit5;  c3 += hit3
    return (c5 / draws, c3 / draws)


def main():
    if len(sys.argv) < 4:
        print(__doc__);  sys.exit(2)
    ex, guard, log = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    here    = os.path.dirname(os.path.abspath(log))
    sums    = [os.path.join(here, os.path.basename(s)) for s in summaryFiles(ex)]
    printed = printedNumbers([log] + sums + sys.argv[4:])
    quoted, scope = quotedNumbers(ex, guard)

    print(f"\nexample      {ex}  (guard at line {guard})")
    print(f"log          {log}")
    print(f"summary      {', '.join(os.path.basename(s) for s in sums) if sums else '(none opened)'}")
    print(f"checked      {scope}")
    print(f"the run printed {len(printed)} distinct numbers;  {len(quoted)} of 4+ digits are checked")
    if not quoted:
        print("\n>> this branch quotes no number of 4+ digits in its own comments, so there is NOTHING TO CHECK.")
        print(">> That is the usual case -- on 28-Aug-2026 only 4 of 38 branches quoted a comparable number --")
        print(">> and it is NOT a finding.")
        return

    meanCand = sum(len(c) for _, c in quoted) / len(quoted)
    r5, r3   = falseMatchRate(printed, meanCand)
    print(f"\n{'quoted':>18}  {'read as':>16}  {'closest printed':>18}  {'rel.diff':>10}  {'abs.diff':>11}  verdict")
    print('-' * 105)
    nRep, nLast, notFound = 0, 0, []
    for s, cand in quoted:
        best, rel, used = None, None, cand[0]
        for v in cand:
            b, r = closest(v, printed)
            if r is not None and (rel is None or r < rel):
                best, rel, used = b, r, v
        if   rel is not None and rel < REPRODUCED:   verdict = 'REPRODUCED';   nRep  += 1
        elif rel is not None and rel < LASTDIGITS:   verdict = 'last digits';  nLast += 1
        else:                                        verdict = 'NOT FOUND';    notFound.append(s)
        absd = abs(abs(used) - best) if best is not None else None
        print(f"{s:>18}  {used:>16.6g}  {(('%.6g' % best) if best is not None else '-'):>18}  "
              f"{('%.2e' % rel) if rel is not None else '-':>10}  "
              f"{('%.3e' % absd) if absd is not None else '-':>11}  {verdict}")
    print('-' * 105)
    print(f">> REPRODUCED to better than {REPRODUCED:g}:  {nRep} of {len(quoted)}."
          f"   Expected by chance here: {r5 * len(quoted):.1f}"
          f"  (false-match rate {r5:.1%} at {meanCand:.1f} candidate scalings each).")
    if nLast:
        print(f">> Agreeing only to {LASTDIGITS:g}:  {nLast}."
              f"   Expected by chance: {r3 * len(quoted):.1f}  (rate {r3:.1%}) -- too weak to count as either.")
    if notFound:
        print(f">> NOT FOUND among the printed numbers:  {len(notFound)}  --  {', '.join(notFound[:12])}"
              f"{' ...' if len(notFound) > 12 else ''}")
        print(">> A NOT FOUND IS A QUESTION, NOT A FINDING, and the tool cannot make it one: a number that moved")
        print(">> and a number that was never a printed result -- a literature value, an input, a ratio, a")
        print(">> threshold -- are indistinguishable from here.  READ THE COMMENT AROUND EACH before filing it.")


main()
