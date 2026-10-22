# Compares the MR_3 electron-impact cross sections of O 2p^4 3P -> 2p^3 3s 3S and 2p^3 3d 3D (initial 3P_J averaged with 2J+1, final J
# summed) for the two exchange models with the digitized experiments. Reads the output of tools/run-ie-oxygen-mr3.jl and the
# *-clean.csv files of tools/digitize-oxygen-*.py from <dir>; writes ics-mr3.png and dcs-3d3D-mr3.png there.
#
#   python3 -I tools/plot-ie-oxygen-mr3.py <dir, e.g. work/o-mr3>

import sys, os, csv, glob, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "work/o-mr3"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
MODELS = [("FurnessMcCarthy", "MR_3, Furness-McCarthy exchange", "#eb6834", "--"),
          ("nonlocal",        "MR_3, non-local exchange",        "#1baf7a", "-")]
EXPT = {"Kanik2001":           ("Kanik et al 2001",           "#2a78d6", "o"),
        "VaughanDoering":      ("Vaughan & Doering 1987/88",  "#7a4fc4", "s"),
        "VaughanDoering8687":  ("Vaughan & Doering 1986/87",  "#7a4fc4", "s"),
        "DoeringYang2001":     ("Doering & Yang 2001",        "#6b6a64", "^"),
        "Johnson2003":         ("Johnson et al 2003",         "#6b6a64", "v"),
        "GulcicekDoering1988": ("Gulcicek & Doering 1988",    "#6b6a64", "D")}

def rows(path):
    with open(path) as f:
        return [r for r in csv.reader(f) if r and not r[0].startswith("#")]

def style(ax, xlabel, ylabel, log=False):
    if log:  ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, linewidth=0.6);  ax.set_axisbelow(True)
    for s in ("top", "right"):    ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):  ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_xlabel(xlabel, color=INK, fontsize=9);  ax.set_ylabel(ylabel, color=INK, fontsize=9)

def expt(ax, path):
    for r in rows(path)[1:]:
        name, x, y, lo, hi = r[0], float(r[1]), float(r[2]), float(r[3]), float(r[4])
        lab, col, mk = EXPT[name]
        yerr = None if math.isnan(lo) else [[y - lo], [hi - y]]
        ax.errorbar([x], [y], yerr=yerr, fmt=mk, color=col, ms=6, capsize=2, elinewidth=1, mec="white", mew=0.8,
                    label=lab if lab not in ax.get_legend_handles_labels()[1] else None, zorder=3)

# Integrated cross sections
ics = {}
for m, *_ in MODELS:
    p = os.path.join(D, "mr3-%s-ics.csv" % m)
    if os.path.exists(p):
        rr = rows(p);  head = [h.strip().split()[0] for h in rr[0]]
        ics[m] = {t: sorted((float(r[1]), float(r[head.index(t)])) for r in rr[1:]) for t in ("3s3S", "3d3D")}
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, (term, title) in zip(axes, [("3s3S", "2p$^4$ $^3$P $\\rightarrow$ 2p$^3$3s $^3$S$^o$"), ("3d3D", "2p$^4$ $^3$P $\\rightarrow$ 2p$^3$3d $^3$D$^o$")]):
    for m, lab, col, ls in MODELS:
        if m in ics and ics[m][term]:
            x, y = zip(*ics[m][term]);  ax.plot(x, y, ls, color=col, lw=2, marker=".", ms=5, label=lab, zorder=2)
    expt(ax, os.path.join(D, "expt", "ics-%s-clean.csv" % term))
    ax.set_xlim(0, 210);  ax.set_ylim(bottom=0)
    style(ax, "electron energy  [eV]", "integrated cross section  [10$^{-18}$ cm$^2$]")
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=7.5)
fig.tight_layout();  fig.savefig(os.path.join(D, "ics-mr3.png"), dpi=150, facecolor="white");  print("wrote ics-mr3.png")

# Differential cross sections, 3d 3D at 30, 50, 100 eV
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=False)
for ax, E in zip(axes, (30, 50, 100)):
    for m, lab, col, ls in MODELS:
        p = os.path.join(D, "mr3-%s-dcs-%s-%geV.csv" % (m, m, E))
        if os.path.exists(p):
            rr = rows(p);  head = [h.strip().split()[0] for h in rr[0]];  k = head.index("3d3D")
            ax.plot([float(r[0]) for r in rr[1:]], [max(float(r[k]), 1e-30) for r in rr[1:]], ls, color=col, lw=2, label=lab, zorder=2)
    expt(ax, os.path.join(D, "expt", "dcs-3d3D-%deV-clean.csv" % E))
    ax.set_xlim(0, 180);  ax.set_xticks(range(0, 181, 30))
    style(ax, "scattering angle  [deg]", "DCS  [10$^{-18}$ cm$^2$/sr]", log=True)
    ax.set_title("3d $^3$D$^o$,  %d eV" % E, fontsize=11, color=INK, loc="left")
axes[0].legend(frameon=False, fontsize=7.5)
fig.tight_layout();  fig.savefig(os.path.join(D, "dcs-3d3D-mr3.png"), dpi=150, facecolor="white");  print("wrote dcs-3d3D-mr3.png")
