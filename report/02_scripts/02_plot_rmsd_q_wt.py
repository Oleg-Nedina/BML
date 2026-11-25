
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[2]
base = root / "02_wt_baseline"
out = root / "report" / "02_utilities"
out.mkdir(parents=True, exist_ok=True)

# 1000 MD steps * 2 fs = 2 ps = 0.002 ns
dt_ns = 0.002

data = {}

for temp in ["300K", "400K"]:
    rmsd = np.loadtxt(base / f"rmsd_{temp}.dat")
    q = np.loadtxt(base / f"q_time_{temp}.dat")

    # Frame 0 is the first DCD frame (t = 0 ns, where Q = 1.0 and RMSD = 0.0)
    time_ns = rmsd[:, 0] * dt_ns

    data[temp] = {
        "time": time_ns,
        "rmsd": rmsd[:, 1],
        "q": q[:, 1],  # Q_all (starts from 1.0 at t=0)
    }

def moving_average(values, window=25):
    # 25 frames correspond to approximately 50 ps
    if len(values) < window:
        return values
    kernel = np.ones(window) / window
    smoothed = np.convolve(values, kernel, mode="valid")
    left = (window - 1) // 2
    right = len(values) - len(smoothed) - left
    return np.pad(smoothed, (left, right), mode="edge")

colors = {"300K": "tab:blue", "400K": "tab:red"}

for metric, ylabel, filename in [
    ("rmsd", "Backbone RMSD (Å)", "WT_RMSD_time.png"),
    ("q", "Fraction of native contacts, Q", "WT_Q_time.png"),
]:
    fig, ax = plt.subplots(figsize=(10, 5))

    for temp in ["300K", "400K"]:
        t = data[temp]["time"]
        y = data[temp][metric]

        ax.plot(t, y, color=colors[temp], alpha=0.2)
        ax.plot(
            t,
            moving_average(y),
            color=colors[temp],
            linewidth=2,
            label=f"WT {temp}",
        )

    ax.axvline(
        1.5,
        linestyle="--",
        color="gray",
        alpha=0.5,
        label="End of 300 K simulation",
    )

    ax.set_xlabel("MD time (ns, approximate)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"WT {metric.upper()} over time")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(out / filename, dpi=250)
    plt.close(fig)

print("Plots saved in:", out)
