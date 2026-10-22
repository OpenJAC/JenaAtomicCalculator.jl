# Compares the MR_3 runs of tools/run-ie-oxygen-mr3.jl with each other and with the digitized experiments: O 2p^4 3P -> 2p^3 3s 3S and
# 2p^3 3d 3D, initial 3P_J averaged with 2J+1, final J summed. Writes into <dir>:
#   compare-ics.png          integrated cross sections, both terms, all runs + experiments
#   compare-dcs-3d3D.png     3d 3D DCS at 30, 50, 100 eV, all runs + experiments
#   compare-dcs-3s3S.png     3s 3S DCS at 30, 50, 100 eV, all runs
#   compare-ics.csv          the integrated cross sections side by side, with the ratio to the first run
#
#   python3 -I tools/plot-ie-oxygen-mr3-compare.py <dir, e.g. O_I/o-mr3>

import sys, os, csv, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "O_I/o-mr3"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
# key, label, ics file, dcs file pattern (energy), colour, line style
RUNS = [("R1-FM", "R1  FM, kappa 20",                             "mr3-FurnessMcCarthy-ics.csv", "mr3-FurnessMcCarthy-dcs-FurnessMcCarthy-%geV.csv", "#9b9a94", "--"),
        ("R1-NL", "R1  non-local, kappa 20",                      "mr3-nonlocal-ics.csv",        "mr3-nonlocal-dcs-nonlocal-%geV.csv",               "#9b9a94", ":"),
        ("R2",    "R2  FM + top-up + threshold shift",            "r2/r2-ics.csv",               "r2/r2-dcs-FurnessMcCarthy-%geV.csv",               "#eb6834", "-"),
        ("R3",    "R3  as R2, kappa 40",                          "r3/r3-ics.csv",               "r3/r3-dcs-FurnessMcCarthy-%geV.csv",               "#2a78d6", "-"),
        ("R4",    "R4  non-local + top-up + threshold shift",     "r4/r4-ics.csv",               "r4/r4-dcs-nonlocal-%geV.csv",                      "#1baf7a", "-")]
EXPT = {"Kanik2001":           ("Kanik et al 2001",           "o", INK,   INK),
        "VaughanDoering":      ("Vaughan & Doering 1987/88",  "s", INK,   "white"),
        "VaughanDoering8687":  ("Vaughan & Doering 1986/87",  "s", INK,   "white"),
        "DoeringYang2001":     ("Doering & Yang 2001",        "^", MUTED, MUTED),
        "Johnson2003":         ("Johnson et al 2003",         "v", MUTED, MUTED),
        "GulcicekDoering1988": ("Gulcicek & Doering 1988",    "D", MUTED, "white")}
TERMS = [("3s3S", "2p$^4$ $^3$P $\\rightarrow$ 2p$^3$3s $^3$S$^o$"), ("3d3D", "2p$^4$ $^3$P $\\rightarrow$ 2p$^3$3d $^3$D$^o$")]

def rows(path):
    with open(path) as f:
        return [r for r in csv.reader(f) if r]

def column(path, term):
    rr = rows(path)
    if len(rr) < 2:  return [], None                                # a run still in progress has only opened its files
    head = [h.strip().split()[0] for h in rr[0]];  k = head.index(term)
    return rr, k

def style(ax, xlabel, ylabel, log=False):
    if log:  ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, linewidth=0.6);  ax.set_axisbelow(True)
    for s in ("top", "right"):    ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):  ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_xlabel(xlabel, color=INK, fontsize=9);  ax.set_ylabel(ylabel, color=INK, fontsize=9)

def expt(ax, path):
    if not os.path.exists(path):  return
    for r in rows(path)[1:]:
        name, x, y, lo, hi = r[0], float(r[1]), float(r[2]), float(r[3]), float(r[4])
        lab, mk, ec, fc = EXPT[name]
        yerr = None if math.isnan(lo) else [[y - lo], [hi - y]]
        ax.errorbar([x], [y], yerr=yerr, fmt=mk, color=ec, mfc=fc, mec=ec, ms=5.5, capsize=2, elinewidth=0.9,
                    label=lab if lab not in ax.get_legend_handles_labels()[1] else None, zorder=4)

def ics(run, term):
    p = os.path.join(D, run[2])
    if not os.path.exists(p):  return []
    rr, k = column(p, term)
    if k is None:  return []
    return sorted((float(r[1]), float(r[k])) for r in rr[1:])

# Integrated cross sections
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
for ax, (term, title) in zip(axes, TERMS):
    for run in RUNS:
        pts = ics(run, term)
        if pts:
            x, y = zip(*pts);  ax.plot(x, y, run[5], color=run[4], lw=1.8, marker=".", ms=5, label=run[1], zorder=3)
    expt(ax, os.path.join(D, "expt", "ics-%s-clean.csv" % term))
    ax.set_xlim(0, 210);  ax.set_ylim(bottom=0)
    style(ax, "electron energy  [eV]", "integrated cross section  [10$^{-18}$ cm$^2$]")
    ax.set_title(title, fontsize=11, color=INK, loc="left");  ax.legend(frameon=False, fontsize=7)
fig.tight_layout();  fig.savefig(os.path.join(D, "compare-ics.png"), dpi=150, facecolor="white")

# Differential cross sections
for term, title in TERMS:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, E in zip(axes, (30, 50, 100)):
        for run in RUNS:
            p = os.path.join(D, run[3] % E)
            if os.path.exists(p):
                rr, k = column(p, term)
                if k is None:  continue
                ax.plot([float(r[0]) for r in rr[1:]], [max(float(r[k]), 1e-30) for r in rr[1:]], run[5], color=run[4], lw=1.8,
                        label=run[1], zorder=3)
        expt(ax, os.path.join(D, "expt", "dcs-%s-%deV-clean.csv" % (term, E)))
        ax.set_xlim(0, 180);  ax.set_xticks(range(0, 181, 30))
        style(ax, "scattering angle  [deg]", "DCS  [10$^{-18}$ cm$^2$/sr]", log=True)
        ax.set_title(title + ",  %d eV" % E, fontsize=10, color=INK, loc="left")
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout();  fig.savefig(os.path.join(D, "compare-dcs-%s.png" % term), dpi=150, facecolor="white")

# Table
with open(os.path.join(D, "compare-ics.csv"), "w") as f:
    f.write("term,energy_eV," + ",".join(r[0] for r in RUNS) + "," + ",".join(r[0] + "/R1-FM" for r in RUNS[1:]) + "   # 1e-18 cm^2\n")
    for term, _ in TERMS:
        data = {r[0]: dict(ics(r, term)) for r in RUNS}
        for E in sorted(set().union(*[d.keys() for d in data.values()])):
            v = [data[r[0]].get(E, float("nan")) for r in RUNS]
            ratio = [x / v[0] if v[0] == v[0] and v[0] else float("nan") for x in v[1:]]
            f.write("%s,%g,%s,%s\n" % (term, E, ",".join("%.4g" % x for x in v), ",".join("%.4f" % x for x in ratio)))
print("wrote compare-ics.png, compare-dcs-3d3D.png, compare-dcs-3s3S.png, compare-ics.csv in", D)
