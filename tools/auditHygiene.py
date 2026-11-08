#!/usr/bin/env python3
"""Measure the Rule 15 code-hygiene criteria over src/.  READ-ONLY -- it changes nothing.

    usage:   python3 tools/auditHygiene.py [--json out.json] [--files]

WHY THIS FILE EXISTS AT ALL.  Five separate hygiene detectors were written against this
codebase and FIVE OF THEM LIED, each by a factor between three and fifty, and each in a way
that looked entirely plausible until a sample was read by eye.  A hygiene count that nobody
can re-derive is worth nothing, so the measurement lives here rather than in a session.
Every calibration below was paid for with a wrong answer; do not simplify one away without
reading a sample of what it changes.

    1. AN ARGUMENT LIST MUST BE SPLIT AT BRACE DEPTH ZERO.  Splitting on every comma reads
       `Array{Tuple{Int64,Float64},1}` as two arguments with no type, and reported 30 % of
       all functions as untyped where the true figure is 1.3 %.
    2. A DOCSTRING HEADER IS COMPARED ON ITS FIRST BACKTICKED NAME, not its last.  Taking
       the last reported 14 header mismatches, every one of them an artifact of a signature
       wrapped over two lines.
    3. A DOCSTRING'S CLOSING TRIPLE QUOTE IS NOT A NEW DEFINITION.  Counting it as one
       reported 4159 definitions with too few blank lines above them; the true figure is 291.
    4. A CONTINUATION LINE IS NOT A BAD INDENT.  A line aligned under an open bracket, or
       under a `using`/`export` list, or a wrapped line of comment prose, is correct where it
       is.  Counting them reported 3282 indent violations where there are, in practice, none.
    5. `Module.` IN A DOCSTRING SIGNATURE IS THE CONVENTION, NOT A DEFECT.  JAC writes
       ``Hfs.Settings(set::Hfs.Settings; ...)`` in the docstring above `function Settings(...)`.
       2017 docstrings differ from their function line for that reason alone.
    6. THE TWO `RacahAlgebra-inc-` FILES ARE OUT OF SCOPE ENTIRELY -- the maintainer's
       decision of 28-Aug-2026.  Their `#` lines are structural to how the sum rules are
       written out.
    7. A `struct` BLOCK MAKES A NESTED `function` AN INNER CONSTRUCTOR, which is correctly
       indented by four.  Counting those reported 10 top-level definitions off column 1 where
       there are 8.
    8. THE BLANK-LINE RULE IS ABOUT THE GAP BETWEEN TWO CODE BLOCKS, so the line above the gap
       must be an `end`.  Counting a gap after a `using`, `include` or `export` line as well
       reported 101 defects that the rule does not cover, and made the sweep look as though it
       had left work behind when it had correctly skipped them.
    9. `=== nothing` AND `!== nothing` ARE THE IDIOMATIC FORM, and are what `isnothing`
       compiles to.  A regex for `== nothing` matches the tail of both and reported 19 guards
       to repair where `src/` holds none at all.
   10. A COMMENT INSIDE A `#= ... =#` BLOCK is part of a deliberately disabled block, which
       Rule 15 keeps verbatim;  it is counted separately rather than offered as work.
   11. A NESTED LOCAL HELPER NEEDS NO DOCSTRING OF ITS OWN.  `module-RadialIntegrals.jl`
       defines the Slater kernel as a two-line `ul(r, s)` inside each parent whose docstring
       already states the integral;  all ten of that file's reported missing docstrings are of
       that kind.  Counting them put the tree-wide figure at 154 against a true 36.
   12. A DOCSTRING MAY WRAP ITS OWN RETURN STATEMENT ACROSS A LINE BREAK, so the text must be
       flattened before the search.  83 do, and counting line by line put the missing-return
       figure at 602 where it is 519 -- and would have appended a SECOND return sentence to
       docstrings that already carried one.
   13. A DOCSTRING'S SIGNATURE MAY SPAN MANY LINES, so the search window is the whole docstring.
       A four-line window saw an opening backtick with no closing one and reported 148 docstrings
       as carrying no signature where 88 do -- and the 60 extra were the very docstrings this
       project had just corrected to reproduce the function's full typed header.
   14. A GENERIC `function name end` DECLARATION HAS NO BODY and therefore no return value and no
       signature of its own;  `module-Basics.jl` holds 103 of them, some with a docstring.  They
       were excluded from the missing-docstring count and not from the others.
   15. THE MAINTAINER'S RULINGS OF 08-Nov-2026 CHANGED WHAT IS A DEFECT HERE, and two counters moved:
       `textBanner` is RETIRED, because every `#####` separator is now kept untouched whatever it carries --
       the old test flagged three good section labels;  and `noReturnStated` is a JUDGEMENT rather than a
       defect, since a return value may be left unsaid where the method's NAME settles it, so there is no
       count to drive to zero.  An appended `##` is now a defect like a standalone one, with the single
       exception in note 16.
   16. AN APPENDED `##` IN A DATA TABLE IS A ROW INDEX.  `module-LSjj-inc.jl` alone carries
       20541 of them, numbering the jj-LS coefficient rows against the published tables.
       They are reported separately from ordinary appended comments and are never swept.
"""
import os, re, sys, json, collections

SRC  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
SKIP = {"module-RacahAlgebra-inc-sumrules.jl", "module-RacahAlgebra-inc-special.jl"}
DATA = {"module-LSjj-inc.jl", "module-LSjj.jl", "module-PeriodicTable.jl",
        "module-Empirical-inc-freefree-gaunt-data.jl", "module-SpinAngularTables-inc-data.jl"}

FUNANY = re.compile(r'^(\s*)function\s+([A-Za-z_][\w.!]*)\s*\(')
FUNCOL = re.compile(r'^function\s+(.+?)\s*$')
GENDEC = re.compile(r'^\s*function\s+[A-Za-z_][\w.!]*\s*end\s*$')
BLOCK  = re.compile(r'^\s*(module|function|for|while|if|let|begin|try|quote|struct|mutable\s+struct|macro)\b')
TRIVIAL = re.compile(r'^(Base\.(show|string|isless|unique|hash|length|getindex|==))')
RETURNED = re.compile(r'is returned|are returned|returns|nothing is returned', re.I)


def docMask(L):
    """-> (inside-docstring flag per line, is-an-opening-delimiter flag per line)."""
    inDoc, doc, opener = False, [False]*len(L), [False]*len(L)
    for i, l in enumerate(L):
        q = l.count('"""')
        if not inDoc and q:
            doc[i] = True
            if q == 1:  inDoc = True;  opener[i] = True
        elif inDoc:
            doc[i] = True
            if q:  inDoc = False
    return doc, opener


def splitArgs(s):
    """Split an argument list at brace, bracket and paren depth zero."""
    out, cur, d = [], "", 0
    for ch in s:
        if   ch in "({[":  d += 1
        elif ch in ")}]":  d -= 1
        if ch in ",;" and d == 0:  out.append(cur);  cur = ""
        else:                      cur += ch
    if cur.strip():  out.append(cur)
    return [a for a in out if a.strip()]


def signature(L, k):
    """The full signature at line k, continued over wrapped lines and CUT at the close of the
    argument list.  A one-line definition puts its body on the same line, so returning the whole
    line reported `determinePolarizationLambda(kind::LeftCircular)   return(  1 )    end` as
    disagreeing with a docstring that states the signature correctly:  7 of 115."""
    sig, t = FUNCOL.match(L[k]).group(1) if FUNCOL.match(L[k]) else L[k].strip(), k
    while sig.count("(") > sig.count(")") and t+1 < len(L):
        t += 1;  sig += " " + L[t].strip()
    q = argParen(sig)
    if q >= 0:
        d = 0
        for i in range(q, len(sig)):
            if   sig[i] == "(":  d += 1
            elif sig[i] == ")":
                d -= 1
                if d == 0:  return sig[:i+1] + whereOf(sig[i+1:])
    return sig


def whereOf(tail):
    """The ` where {...}` clause of a signature, which belongs to it, or an empty string."""
    m = re.match(r'\s*(where\s+\S+)', tail)
    return " " + m.group(1) if m else ""


def argParen(s):
    """Index of the parenthesis that opens the ARGUMENT list, or -1.

    `Base.:(==)(confa::Configuration, confb::Configuration)` carries an earlier parenthesis that
    belongs to the operator's NAME, and splitting at the first one reads the name as `Base.:` and
    the argument list as `(==)(confa...)`.  That was symmetric between docstring and code and so
    hid, until the name was also required to be an identifier:  17 `Base.:(==)` methods were then
    reported as documenting no signature at all.
    """
    for i, ch in enumerate(s):
        if ch == "(" and (i == 0 or s[i-1] != ":"):  return i
    return -1


def splitSig(s):
    """-> (bare name, normalised argument list with `Module.` prefixes dropped), or None."""
    p = argParen(s)
    if p < 0:  return None
    name = s[:p].strip()
    name = name.split(".")[-1] if "." in name else name
    args = re.sub(r'\b[A-Za-z_]\w*\.(?=[A-Z])', '', s[p:].split(" where")[0])
    return name, re.sub(r'\s+', '', args)


def headersOf(head):
    """-> (candidate signatures, headers).  Both are backticked fragments that START a line and name
       a function, joined across the backtick pairs they span;  a HEADER is additionally followed by the
       `...` that opens a description.  The looser set answers "which signature does this docstring
       state", the stricter one "how many things does it document", and they must stay separate:
       requiring the `...` of the first reported 56 docstrings as carrying no signature where 23 do,
       while not requiring it of the second reported 106 families where there are 24."""
    spans = [(m.start(1), m.end(1), m.group(1)) for m in re.finditer(r'`([^`]+)`', head, re.S)]
    cands, heads, i = [], [], 0
    while i < len(spans):
        b, e, txt = spans[i]
        ## THREE TESTS, AND ALL THREE ARE NEEDED.  Requiring only "starts a line" read a tuple, a
        ## formula and a cited call as headers and reported 106 family docstrings where there are 37:
        ## `(orbitals, chemMu, ...)` opening a returned NamedTuple, `n = (2 temp)^(3/2) ...` opening a
        ## displayed equation, and `minimum(g.meanEnergy for g in groupsPerStage[i+1])` opening a
        ## sentence all start a line and all carry a parenthesis.
        bol = head.rfind("\n", 0, b) + 1
        if head[bol:b-1].strip().strip('"') != "" or "(" not in txt:  i += 1;  continue
        ## (2) THE NAME MUST BE A DOTTED IDENTIFIER, which an equation and a tuple are not.
        q = argParen(txt)
        if q < 0 or not re.match(r'^[A-Za-z_]\w*(\.([A-Za-z_]\w*|:\([^)]*\)|:\S+))*[!?]?$', txt[:q].strip()):
            i += 1;  continue
        ## A WRAPPED SIGNATURE IS CONTINUED ALONG THE DOCSTRING'S OWN LINES, backticked or not, and
        ## the walk STOPS at a blank line or at the line that opens the description.  Joining across
        ## backtick SPANS instead, which ignores the text and the blank lines between them, mangled 13
        ## of 122 reported mismatches:  `Hfs.Settings(; calcM1::Bool=true,` closes its backtick in
        ## mid-signature and the remainder of the keyword list is unbackticked, so the walk ended at
        ## the `settings::Hfs.Settings` of the return sentence four lines below.
        acc, endPos = txt, e+1
        while acc.count("(") > acc.count(")"):
            eol  = head.find("\n", endPos)
            if eol < 0:  eol = len(head)
            rest = head[endPos:eol]
            if rest.strip():
                acc += rest.replace("`", "");  endPos = eol
                if acc.count("(") <= acc.count(")"):  break
            if eol >= len(head):  break
            nEnd = head.find("\n", eol+1)
            if nEnd < 0:  nEnd = len(head)
            if not head[eol+1:nEnd].strip() or head[eol+1:nEnd].strip().startswith(("...", "+")):  break
            endPos = eol + 1
        j = i
        while j+1 < len(spans) and spans[j+1][0] < endPos:  j += 1
        ## (3) A HEADER IS FOLLOWED BY THE `...` THAT OPENS ITS DESCRIPTION, on the rest of its own
        ## line or on the next non-blank one.  This is what separates a header from a cited call.
        eol  = head.find("\n", endPos)
        tail = head[endPos : eol if eol >= 0 else len(head)].strip().strip("`").strip()
        if not tail:
            rest = head[eol+1:] if eol >= 0 else ""
            tail = next((l.strip() for l in rest.split("\n") if l.strip()), "")
        cands.append(acc)
        if tail.startswith("..."):  heads.append(acc)
        i = j+1
    return cands, heads


def disabledMask(L):
    """-> flag per line: inside a `#= ... =#` block, which Rule 15 keeps verbatim."""
    out, depth = [False]*len(L), 0
    for i, l in enumerate(L):
        out[i] = depth > 0
        if depth:
            depth = max(0, depth + l.count("#=") - l.count("=#"));  continue
        if "#=" in l and "=#" not in l.split("#=", 1)[1]:  depth = 1
    return out


def audit(path):
    fn = os.path.basename(path)
    L  = open(path, errors="replace").read().split("\n")
    doc, opener = docMask(L)
    dis = disabledMask(L)
    r = collections.Counter()
    r["lines"] = len(L)
    sites = collections.defaultdict(list)

    # ---------------------------------------------------------------- comments and whitespace
    depth = 0
    for i, l in enumerate(L):
        s = l.strip()
        if len(l) > 140:
            if doc[i] or s.startswith("#"):  r["longProse"] += 1;  sites["longProse"].append(i+1)
            elif fn in DATA:                 r["longData"]  += 1
            else:                            r["longCode"]  += 1
        if doc[i]:  continue
        if dis[i]:  r["inDisabledBlock"] += 1;  continue
        if re.match(r'^\s*#\s*$', l):        r["emptyHash"] += 1;  sites["emptyHash"].append(i+1)
        if re.match(r'^\s*#{5,}\s*$', l):    r["hashRule"]  += 1
        elif re.match(r'^\s*#{5,}', l) and re.sub(r'[#\s]', '', l):
            r["textBanner"] += 1;  sites["textBanner"].append(i+1)
        elif re.match(r'^\s*##[^#]', l):     r["standalone2"] += 1;  sites["standalone2"].append(i+1)
        elif re.search(r'\S\s+##[^#]', l) and not s.startswith("#"):
            r["appended2Data" if fn in DATA else "appended2"] += 1
        if re.search(r'(?<![=!])==\s*nothing', l):
            r["eqNothing"] += 1;  sites["eqNothing"].append(i+1)
        if s and not s.startswith("#") and depth == 0:
            ind = len(l) - len(l.lstrip())
            if ind % 4 and ind < 40:  r["indentOdd"] += 1
        if not s.startswith("#"):
            code = re.sub(r'".*?"', '', l)
            depth = max(0, depth + code.count("(") + code.count("[") - code.count(")") - code.count("]"))

    # ---------------------------------------------------------------- functions
    stack = []
    for i, l in enumerate(L):
        if GENDEC.match(l):  r["genericDecl"] += 1;  continue
        m = FUNANY.match(l)
        if m:
            r["funcs"] += 1
            name = m.group(2)
            r["qualified" if "." in name else "bare"] += 1
            nested = bool(m.group(1)) and ("function" in stack or "struct" in stack)
            if not m.group(1):                 pass
            elif nested:                       r["nestedFunc"] += 1
            else:  r["indentedTop"] += 1;      sites["indentedTop"].append(i+1)
            args = splitArgs(signature(L, i).split("(", 1)[-1].rsplit(")", 1)[0])
            if args and not all("::" in a or "=" in a or "..." in a for a in args):
                r["untyped"] += 1;  sites["untyped"].append(i+1)
            # docstring immediately above?
            j = i-1
            while j >= 0 and L[j].strip() == "":  j -= 1
            if j >= 0 and L[j].strip().endswith('"""') and doc[j]:
                pass
            elif TRIVIAL.match(name):  r["noDocTrivial"] += 1
            elif nested:               r["noDocNested"] += 1
            else:                      r["noDoc"] += 1;  sites["noDoc"].append(i+1)
        if doc[i] or l.strip().startswith("#"):  continue
        if BLOCK.match(l) and not re.search(r'\bend\b\s*$', l) and not re.match(r'^\s*\w+\s*\(.*\)\s*=', l):
            stack.append("function" if re.match(r'^\s*function\b', l) else
                         "struct"   if re.match(r'^\s*(mutable\s+)?struct\b', l) else "other")
        elif re.match(r'^\s*end\b', l) and stack:
            stack.pop()

    # ---------------------------------------------------------------- docstring bodies
    i = 0
    while i < len(L):
        if not (opener[i] and L[i].startswith('"""')):  i += 1;  continue
        j = i+1
        while j < len(L) and '"""' not in L[j]:  j += 1
        k = j+1
        while k < len(L) and L[k].strip() == "":  k += 1
        ## A GENERIC `function name end` DECLARATION HAS NO BODY, so it can state no return value
        ## and has no signature of its own.  `module-Basics.jl` holds 103 of them, several carrying
        ## a docstring, and counting them put the missing-return figure one too high per docstring.
        if k < len(L) and FUNCOL.match(L[k]) and not GENDEC.match(L[k]):
            r["funcDocs"] += 1
            ## THE WHITESPACE IS NORMALISED BEFORE THE SEARCH.  A docstring may wrap its own return
            ## statement across a line break -- "Nothing is\n              returned." -- and 83 of
            ## them do;  searching line by line put the figure at 602 where it was 519.
            body = re.sub(r'\s+', ' ', "\n".join(L[i:j+1]))
            if not RETURNED.search(body):  r["noReturnStated"] += 1;  sites["noReturnStated"].append(k+1)
            ## THE WINDOW IS THE WHOLE DOCSTRING, not its first four lines.  A copy-constructor
            ## signature reproduced from the code can run to eight lines, and a four-line window
            ## then sees an opening backtick with no closing one and reports "no signature at
            ## all" -- it put the figure at 148 where it is 88, and every one of the 60 extra was
            ## a docstring this project had just CORRECTED.
            cands, heads = headersOf("\n".join(L[i:j]))
            if not cands:
                r["noSigInDoc"] += 1;  sites["noSigInDoc"].append(k+1)
            else:
                ## A DOCSTRING CARRYING MORE THAN ONE HEADER DOCUMENTS MORE THAN ONE THING, which
                ## Rule 15.3 does not allow:  a reader arriving at one method must find ITS
                ## description attached to it.  Two shapes occur -- several methods sharing one
                ## docstring, and one method whose two call forms are written as two headers --
                ## and both end at one docstring per method whose header matches the code.
                if len(heads) > 1:  r["docFamily"] += 1;  sites["docFamily"].append(k+1)
                b = splitSig(signature(L, k))
                if b and not any(splitSig(h) == b for h in cands  if splitSig(h)):
                    r["sigMismatch"] += 1;  sites["sigMismatch"].append(k+1)
        i = j+1

    # ---------------------------------------------------------------- blank-line separation
    for i, l in enumerate(L):
        if GENDEC.match(l) or (doc[i] and not opener[i]):  continue
        isDoc = opener[i] and l.startswith('"""')
        isFun = bool(FUNCOL.match(l))
        if not (isDoc or isFun):  continue
        if isFun and i and L[i-1].strip().endswith('"""'):  continue
        n, j = 0, i-1
        while j >= 0 and L[j].strip() == "":  n += 1;  j -= 1
        if j < 0 or re.match(r'^\s*#', L[j]) or L[j].strip().startswith("export"):  continue
        if not re.match(r'^\s*end\b', L[j]):  r["sepNotAfterEnd"] += 1;  continue
        if   n == 0:  r["sep0"] += 1;  sites["sep0"].append(i+1)
        elif n == 1:  r["sep1"] += 1;  sites["sep1"].append(i+1)
        else:         r["sep2plus"] += 1
    return r, sites


ROWS = [("funcs",           "function definitions"),
        ("genericDecl",     "generic `function name end` declarations"),
        ("funcDocs",        "docstrings attached to a function"),
        ("", ""),
        ("standalone2",     "standalone `##` comment lines        -> `#`"),
        ("appended2",       "appended `##` on a code line        -> `#`"),
        ("emptyHash",       "empty `#` lines                      -> a blank line"),
        ("sep1",            "definitions with only 1 blank line   -> 2"),
        ("sep0",            "definitions with no blank line above"),
        ("sepNotAfterEnd",  "   gaps not after an `end`, so out of scope"),
        ("longProse",       "docstring/comment lines over 140     -> re-flow"),
        ("eqNothing",       "`== nothing` guards                  -> `isnothing`"),
        ("indentedTop",     "top-level functions not at column 1"),
        ## `textBanner` IS RETIRED, 08-Nov-2026:  the maintainer keeps EVERY `#####` separator untouched,
        ## whatever it carries, so there is nothing here to report.  It is counted below, with the other
        ## things that are kept, purely so the number does not look lost.  What WOULD be worth reporting is
        ## the narrow shape Rule 15.7 now names -- a separator block ENCLOSING A PARAGRAPH of prose -- and
        ## that is a different test from "the line carries text", which flagged three good labels.
        ("inDisabledBlock",  "   comment lines inside `#= ... =#`, kept"),
        ("", ""),
        ("noDoc",           "functions with no docstring"),
        ("noDocTrivial",    "   also: trivial Base.* methods (exempt)"),
        ("noDocNested",     "   also: nested local helpers (exempt)"),
        ("noReturnStated",  "   docstrings not stating a return value (JUDGEMENT, not a defect)"),
        ("sigMismatch",     "docstring signatures disagreeing with the code"),
        ("docFamily",       "docstrings carrying more than one header"),
        ("noSigInDoc",      "docstrings carrying no signature at all"),
        ("untyped",         "functions with an argument lacking a type"),
        ("", ""),
        ("textBanner",      "`#####` separators carrying a label (kept, Rule 15.7)"),
        ("nestedFunc",      "nested helper functions (correctly indented)"),
        ("hashRule",        "pure `#####` separators (kept)"),
        ("longCode",        "code lines over 140 (alignment, not prose)"),
        ("longData",        "data-table lines over 140 (kept)"),
        ("appended2Data",   "appended `##` row indices in a data table (kept, the one exception)"),
        ("indentOdd",       "indents not a multiple of 4 (see note 4)"),
        ("qualified",       "`function Module.name(` definitions"),
        ("bare",            "`function name(` definitions")]


def main():
    want  = "--files" in sys.argv
    jpath = sys.argv[sys.argv.index("--json")+1] if "--json" in sys.argv else None
    files = sorted(f for f in os.listdir(SRC) if f.endswith(".jl") and f not in SKIP)
    per, allSites, tot = {}, {}, collections.Counter()
    for f in files:
        r, s = audit(os.path.join(SRC, f))
        per[f] = dict(r);  allSites[f] = {k: v for k, v in s.items()}
        tot.update(r)
    print("src/ hygiene, %d files, %d lines   (the two RacahAlgebra-inc- files are out of scope)"
          % (len(files), tot["lines"]))
    print("-" * 78)
    for key, label in ROWS:
        if not key:  print();  continue
        n = tot[key]
        nf = sum(1 for f in files if per[f].get(key))
        print("   %-52s %6d  in %3d files" % (label, n, nf))
    if want:
        print("\n" + "-" * 78)
        for key, label in ROWS:
            if not key or not tot[key]:  continue
            top = sorted(((per[f].get(key, 0), f) for f in files), reverse=True)[:6]
            print("   %-46s %s" % (label, ", ".join("%s %d" % (f.replace("module-", "").replace(".jl", ""), n)
                                                    for n, f in top if n)))
    if jpath:
        json.dump({"totals": dict(tot), "perFile": per, "sites": allSites}, open(jpath, "w"), indent=0)
        print("\n   per-file counts and line numbers written to %s" % jpath)


main()
