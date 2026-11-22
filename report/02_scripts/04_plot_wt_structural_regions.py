from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "02_wt_baseline"
OUT = ROOT / "report" / "02_utilities"
OUT.mkdir(parents=True, exist_ok=True)

r300 = np.loadtxt(DATA / "rmsf_300K.dat")
r400 = np.loadtxt(DATA / "rmsf_400K.dat")
q = np.genfromtxt(DATA / "q_per_residue_400K.dat", dtype=None, comments="#", encoding="utf-8")
res = np.array([r[0] for r in q], dtype=int)
name = np.array([r[1] for r in q])
contacts = np.array([r[2] for r in q], dtype=int)
qi = np.array([r[3] for r in q], dtype=float)
masked = np.array([r[4] for r in q], dtype=int)
assert np.array_equal(res, r300[:, 0].astype(int)) and np.array_equal(res, r400[:, 0].astype(int))

regions = [
    ("N-terminal", 1, 27),
    ("Ankyrin repeat 1", 28, 60),
    ("Ankyrin repeat 2", 61, 93),
    ("C-terminal", 94, 118),
]

fig, axes = plt.subplots(2, 1, figsize=(13.5, 8.7), sharex=True)
for ax in axes:
    for k, (label, lo, hi) in enumerate(regions):
        ax.axvspan(lo-0.5, hi+0.5, color=("#dfeaf2" if k % 2 == 0 else "#e9e2f4"), alpha=0.50, zorder=0)
    for sep in (27.5, 60.5, 93.5):
        ax.axvline(sep, color="0.45", ls=":", lw=1)
    ax.grid(axis="y", alpha=0.25)
    ax.set_xlim(0.5, 118.5)

axes[0].plot(res, r300[:,1], lw=1.4, color="#2778ad", label="300 K (1.5 ns)")
axes[0].plot(res, r400[:,1], lw=1.6, color="#c85c28", label="400 K (2.5 ns)")
axes[0].set_ylabel("Cα RMSF (Å)")
axes[0].set_ylim(bottom=0)
axes[0].legend(loc="upper right")
fig.suptitle("Myotrophin WT — fluctuations and native-contact retention", fontsize=14, y=0.985)

axes[1].plot(res, qi, lw=1.5, color="#246c6c", marker=".", ms=3,
             label="Qᵢ, final 30% at 400 K")
axes[1].scatter(res[masked == 1], qi[masked == 1], color="#8a609c", s=32, marker="x",
                label="Masked positions (screening exclusion)", zorder=4)
axes[1].axhline(0.60, color="#9f5050", ls="--", lw=1.1,
                label="Exploratory Qᵢ = 0.60")
axes[1].set_ylim(0.25, 1.0)
axes[1].set_ylabel("Mean retained-contact fraction Qᵢ")
axes[1].set_xlabel("Residue number")
axes[1].legend(loc="lower right", fontsize=9)

for j, (label, lo, hi) in enumerate(regions):
    axes[0].text((lo + hi) / 2, 1.015, label, transform=axes[0].get_xaxis_transform(),
                 ha="center", va="bottom", fontsize=9.5, fontweight="bold")

for site, label in [(37, "K37"), (103, "T103")]:
    i = int(np.flatnonzero(res == site)[0])
    axes[0].scatter([site], [r400[i,1]], color="black", s=32, zorder=6)
    axes[0].annotate(label, (site, r400[i,1]), xytext=(site + 2, r400[i,1] + 0.35),
                     arrowprops={"arrowstyle":"-", "color":"0.3"}, fontsize=9)
    axes[1].scatter([site], [qi[i]], color="black", s=35, zorder=6)

fig.text(.5, .015,
         "RMSF windows differ (300 K: 1.5 ns; 400 K: 2.5 ns). Qᵢ uses only the final 30% of 400 K.",
         ha="center", fontsize=9, color="0.35")
fig.tight_layout(rect=(0, .045, 1, .895))
fig.savefig(OUT / "WT_structural_regions_corrected.png", dpi=240)
plt.close(fig)

rows = []
for label, lo, hi in regions:
    in_reg = (res >= lo) & (res <= hi)
    eligible = in_reg & (masked == 0) & (res >= 4) & (res <= 115) & (contacts >= 30)
    low = eligible & (qi < 0.60)
    vals = {
        "Region": label,
        "Residues": f"{lo}-{hi}",
        "Mean_Qi_400K_final30": round(float(qi[in_reg].mean()), 4),
        "Contact_weighted_Qi_400K_final30": round(float(np.average(qi[in_reg], weights=contacts[in_reg])), 4),
        "RMSF_300K_mean_A": round(float(r300[in_reg, 1].mean()), 4),
        "RMSF_400K_mean_A": round(float(r400[in_reg, 1].mean()), 4),
        "Eligible_residues": int(eligible.sum()),
        "Eligible_Qi_under_0.60": int(low.sum()),
    }
    rows.append(vals)

with (OUT / "WT_structural_region_summary.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(row)
print("Saved to:", OUT)
