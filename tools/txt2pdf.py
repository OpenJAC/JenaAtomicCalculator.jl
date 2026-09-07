#!/usr/bin/env python3
"""Plain ASCII text -> PDF, monospaced, 11 pt, A4, no external packages.

    python3 tools/txt2pdf.py  <input.txt>  <output.pdf>

Lives in tools/ so that every application folder can reach it without a copy being carried from
one apps/ directory into another, which the one-application-at-a-time rule forbids and which
happened once before this file was moved here (7-Sep-2026).

Assumes pure ASCII input: any other byte is replaced rather than raising, so a stray character
gives a visible artefact instead of a failed build.  55 lines fit on a page at these margins, so
a page count can be predicted from the line count before rendering.
"""
import sys, zlib

FONT, SIZE, LEAD = "Courier", 11.0, 13.2
PAGE_W, PAGE_H   = 595.28, 841.89          # A4 in points
LEFT, TOP, BOT   = 48.0, 56.0, 48.0
LINES = int((PAGE_H - TOP - BOT) / LEAD)   # lines that fit on one page

def esc(s):
    return s.replace("\\", r"\\").replace("(", r"\(").replace(")", r"\)")

def main(src, dst):
    # rstrip: a trailing newline would otherwise create one extra, empty page.
    raw = open(src, encoding="ascii", errors="replace").read().rstrip("\n").split("\n")
    text = [ln.replace("\t", "    ") for ln in raw]
    pages = [text[i:i+LINES] for i in range(0, len(text), LINES)] or [[""]]

    objs, streams = [], []
    for pg in pages:
        out = ["BT", "/F1 %.1f Tf" % SIZE, "%.2f TL" % LEAD,
               "1 0 0 1 %.2f %.2f Tm" % (LEFT, PAGE_H - TOP)]
        for ln in pg:
            out.append("(%s) Tj T*" % esc(ln))
        out.append("ET")
        streams.append(zlib.compress(("\n".join(out)).encode("latin-1", "replace")))

    n_pg = len(pages)
    # 1 catalog, 2 pages tree, 3 font, then per page: page object and content stream
    objs.append("<< /Type /Catalog /Pages 2 0 R >>")
    kids = " ".join("%d 0 R" % (4 + 2*i) for i in range(n_pg))
    objs.append("<< /Type /Pages /Count %d /Kids [%s] >>" % (n_pg, kids))
    objs.append("<< /Type /Font /Subtype /Type1 /BaseFont /%s >>" % FONT)
    for i in range(n_pg):
        objs.append("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] "
                    "/Resources << /Font << /F1 3 0 R >> >> /Contents %d 0 R >>"
                    % (PAGE_W, PAGE_H, 5 + 2*i))
        objs.append(("<< /Length %d /Filter /FlateDecode >>" % len(streams[i]), streams[i]))

    buf, offsets = bytearray(b"%PDF-1.4\n"), []
    for k, ob in enumerate(objs, start=1):
        offsets.append(len(buf))
        if isinstance(ob, tuple):
            buf += ("%d 0 obj\n%s\nstream\n" % (k, ob[0])).encode()
            buf += ob[1] + b"\nendstream\nendobj\n"
        else:
            buf += ("%d 0 obj\n%s\nendobj\n" % (k, ob)).encode()
    start = len(buf)
    buf += ("xref\n0 %d\n" % (len(objs)+1)).encode() + b"0000000000 65535 f \n"
    for off in offsets:
        buf += ("%010d 00000 n \n" % off).encode()
    buf += ("trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n"
            % (len(objs)+1, start)).encode()
    open(dst, "wb").write(bytes(buf))
    print("%s -> %s :  %d lines, %d pages, %.1f pt %s"
          % (src, dst, len(text), n_pg, SIZE, FONT))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
