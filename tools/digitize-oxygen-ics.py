# Digitizes experimental points (marker centre and error-bar ends) from a matplotlib PNG with LINEAR axes, as the ICS plots in
# O_I/AS1_S3/plots and O_I/oxygen/plots. Both scales come from the major tick marks outside the spines, whose values are given as
# first value and step; the points are found by marker colour. Error-bar ends are approximate where markers of two sets overlap.
#
#   python3 -I tools/digitize-oxygen-ics.py <png> <x0value> <xstep> <y0value> <ystep> <out.csv> name=R,G,B [name=R,G,B ...]

import sys
import numpy as np
from PIL import Image

png = sys.argv[1];   xv0, xst, yv0, yst = map(float, sys.argv[2:6]);   out = sys.argv[6]
SETS = {a.split("=")[0]: tuple(int(v) for v in a.split("=")[1].split(",")) for a in sys.argv[7:]}
img = np.asarray(Image.open(png).convert("RGB")).astype(int)
H, W, _ = img.shape
dark = img.sum(axis=2) < 150

colsum = dark.sum(axis=0);   rowsum = dark.sum(axis=1)
x0 = int(np.argmax(colsum[: W // 2]));   x1 = W // 2 + int(np.argmax(colsum[W // 2:]))
y1 = H // 2 + int(np.argmax(rowsum[H // 2:]));   y0 = int(np.argmax(rowsum[: H // 2]))

# Major ticks: dark marks 4 pixels outside the left and the bottom spine.
yt = sorted([y for y in range(y0, y1 + 1) if dark[y, x0 - 4]], reverse=True)
xt = sorted([x for x in range(x0, x1 + 1) if dark[y1 + 4, x]])
merge = lambda v: [v[i] for i in range(len(v)) if i == 0 or abs(v[i] - v[i - 1]) > 1]
yt, xt = merge(yt), merge(xt)
yval = lambda y: yv0 + yst * (yt[0] - y) / float(np.mean(-np.diff(yt)))
xval = lambda x: xv0 + xst * (x - xt[0]) / float(np.mean(np.diff(xt)))

res = []
for name, c in SETS.items():
    m = (np.abs(img - np.array(c)).sum(axis=2) < 60)
    m[:, : x0 + 1] = False;  m[:, x1:] = False;  m[: y0 + 1, :] = False;  m[y1:, :] = False
    core = m.copy()
    for _ in range(2):   # erode twice: error-bar lines vanish, marker discs survive
        core = core & np.roll(core, 1, 0) & np.roll(core, -1, 0) & np.roll(core, 1, 1) & np.roll(core, -1, 1)
    seen = np.zeros_like(core)
    for y, x in zip(*np.nonzero(core)):
        if seen[y, x]:  continue
        stack, pts = [(y, x)], []
        seen[y, x] = True
        while stack:
            a, b = stack.pop();  pts.append((a, b))
            for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                p, q = a + da, b + db
                if core[p, q] and not seen[p, q]:  seen[p, q] = True;  stack.append((p, q))
        cy = np.mean([p[0] for p in pts]);  cx = np.mean([p[1] for p in pts])
        if cx > 0.55 * W and cy < 0.45 * H:  continue                 # the legend's own markers (upper right)
        col = int(round(cx));  top = int(cy);  bot = int(cy)
        while top > y0 + 1 and m[top - 1, col - 1: col + 2].any():  top -= 1
        while bot < y1 - 1 and m[bot + 1, col - 1: col + 2].any():  bot += 1
        res.append((name, xval(cx), yval(cy), yval(bot), yval(top), cx, cy))

with open(out, "w") as f:
    f.write("set,energy_eV,ics_1e-18cm2,low,high,px,py\n")
    for r in sorted(res, key=lambda r: (r[0], r[1])):  f.write("%s,%.2f,%.4g,%.4g,%.4g,%.1f,%.1f\n" % r)
print("frame x %d-%d y %d-%d; x ticks %s; y ticks %s" % (x0, x1, y0, y1, xt, yt))
for r in sorted(res, key=lambda r: (r[0], r[1])):  print("%-14s %7.2f  %7.3f  [%7.3f, %7.3f]  px (%.0f,%.0f)" % r)
