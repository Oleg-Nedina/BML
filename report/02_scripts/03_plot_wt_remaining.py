from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "02_wt_baseline"
OUT = ROOT / "report" / "02_utilities"
OUT.mkdir(parents=True, exist_ok=True)
DT_NS = 0.002  # 2 ps between saved frames: verify using DCD/log if available


def md_series(variable, temp, column=1):
    """Discard reference PDB (frame 0); return approximate MD time and values."""
    a = np.loadtxt(DATA / f"{variable}_{temp}.dat")
    a = a[a[:, 0] >= 1]
    return (a[:, 0] - 1) * DT_NS, a[:, column]


def add_series(ax, variable, temp, column=1, multiplier=1.0):
    t, values = md_series(variable, temp, column)
    values = values * multiplier
    (line,) = ax.plot(t, values, linewidth=0.8, alpha=0.35, label=f"{temp} raw")
    # 25 data points = ~50 ps, purely for visualization
    window = 25
    if len(values) >= window:
        smooth = np.convolve(values, np.ones(window) / window, mode="valid")
        offset = window // 2
        ax.plot(t[offset:-offset], smooth, color=line.get_color(),
                linewidth=2, label=f"{temp} (50 ps moving mean)")
    else:
        line.set_alpha(1.0)
    return t, values


def finish_time_axes(ax):
    ax.set_xlabel("Approximate MD time (ns)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncol=2)


# FIGURE 1: global expansion and solvent accessibility
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
for temp in ("300K", "400K"):
    add_series(axes[0], "rgyr", temp)
    add_series(axes[1], "sasa", temp)
axes[0].set_ylabel("Radius of gyration (Å)")
axes[0].set_title("WT compactness")
axes[1].set_ylabel("SASA (Å²)")
axes[1].set_title("WT solvent exposure")
for ax in axes:
    finish_time_axes(ax)
fig.tight_layout()
fig.savefig(OUT / "WT_Rg_SASA.png", dpi=250)
plt.close(fig)

# FIGURE 2: secondary structure and hydrogen-bond counts
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
for temp in ("300K", "400K"):
    add_series(axes[0], "alpha", temp, multiplier=100)
    add_series(axes[1], "hbonds", temp)
axes[0].set_ylabel("Residues assigned alpha helix (%)")
axes[0].set_title("Alpha-helical content (VMD secondary-structure assignment)")
axes[1].set_ylabel("Detected intraprotein H-bonds per frame")
axes[1].set_title("Hydrogen-bond count (VMD Hbonds plugin settings)")
for ax in axes:
    finish_time_axes(ax)
fig.tight_layout()
fig.savefig(OUT / "WT_Alpha_Hbonds.png", dpi=250)
plt.close(fig)

# FIGURE 3: script-defined regions and their contacts with the full protein
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, sharey=True)
region_names = [(4, "1–34"), (5, "35–67"), (6, "68–100"), (7, "101–118")]
for ax, temp in zip(axes, ("300K", "400K")):
    a = np.loadtxt(DATA / f"q_time_{temp}.dat")
    a = a[a[:, 0] >= 1]
    t = (a[:, 0] - 1) * DT_NS
    for column, name in region_names:
        ax.plot(t, a[:, column], linewidth=1.0,
                label=f"Region {name}")
    ax.set_title(f"Regional contact retention — {temp}")
    ax.set_ylabel("Fraction of reference contacts")
    ax.set_ylim(0.5, 0.9)
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
axes[-1].set_xlabel("Approximate MD time (ns)")
fig.tight_layout()
fig.savefig(OUT / "WT_Q_regions.png", dpi=250)
plt.close(fig)
