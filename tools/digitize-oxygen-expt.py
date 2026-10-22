# Digitizes the experimental points (marker centre and error-bar ends) from a matplotlib PNG of a DCS on a log y axis with x from 0 to 180
# degrees, as in O_I/AS1_S3/plots/DCS_3d3D_*eV_AS1_S3_vs_MR3.png. The axes frame is found from the dark spines, the y scale from the
# major tick marks on the left spine, and the points by their marker colour.
#
#   python3 -I tools/digitize-oxygen-expt.py <png> <exponent of the lowest major tick, e.g. -2> <out.csv>

import sys, math
import numpy as np
from PIL import Image

png, dec0, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
img = np.asarray(Image.open(png).convert("RGB")).astype(int)
H, W, _ = img.shape
dark = img.sum(axis=2) < 150

# Spines: the longest dark column on the left half and the longest dark row in the lower half.
colsum = dark.sum(axis=0);   rowsum = dark.sum(axis=1)
x0 = int(np.argmax(colsum[: W // 2]));   x1 = W // 2 + int(np.argmax(colsum[W // 2:]))
y1 = H // 2 + int(np.argmax(rowsum[H // 2:]));   y0 = int(np.argmax(rowsum[: H // 2]))

# Y scale from the major ticks: the dark marks that reach 3-5 pixels left of the spine (minor ticks reach only 2).
ticks = [y for y in range(y0, y1 + 1) if dark[y, x0 - 4]]
decades = sorted(ticks, reverse=True)                                # bottom (= 10^dec0) to top
pixPerDec = float(np.mean(-np.diff(decades)))
lines = ticks
ylog = lambda y: dec0 + (decades[0] - y) / pixPerDec
xdeg = lambda x: 180.0 * (x - x0) / (x1 - x0)

SETS = {"Kanik2001": (31, 119, 180), "VaughanDoering": (148, 103, 189)}
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
        if cx > 0.55 * W and cy < 0.35 * H:  continue                 # the legend's own markers (upper right)
        col = int(round(cx))
        top = int(cy);  bot = int(cy)
        while top > y0 + 1 and m[top - 1, col - 1: col + 2].any():  top -= 1
        while bot < y1 - 1 and m[bot + 1, col - 1: col + 2].any():  bot += 1
        res.append((name, xdeg(cx), 10 ** ylog(cy), 10 ** ylog(bot), 10 ** ylog(top), cx, cy))

with open(out, "w") as f:
    f.write("set,theta_deg,dcs_1e-18cm2sr,low,high,px,py\n")
    for r in sorted(res, key=lambda r: (r[0], r[1])):
        f.write("%s,%.2f,%.4g,%.4g,%.4g,%.1f,%.1f\n" % r)
print("frame x %d-%d y %d-%d; major ticks %d; decades at %s; px/decade %.1f" % (x0, x1, y0, y1, len(lines), [round(d) for d in decades], pixPerDec))
for r in sorted(res, key=lambda r: (r[0], r[1])):  print("%-15s %6.1f  %9.4g  [%9.4g, %9.4g]  px (%.0f,%.0f)" % r)
