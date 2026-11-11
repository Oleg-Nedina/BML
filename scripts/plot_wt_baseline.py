#!/usr/bin/env python3
"""
plot_wt_baseline.py
Fase 2 (ROAD2.md - Passo 2.4):
Genera 2 figure complete (300 DPI) per la presentazione in 06_presentation/:
  1. fig5_wt_baseline_plots.png:
     - A: Backbone RMSD(t) 300 K vs 400 K
     - B: Q-value nel tempo (all-all, hydrophilic-all e per Repeat ANK1..4)
     - C: RMSF(CA) per residuo 300 K vs 400 K
     - D: <Q_i> per singolo residuo a 400 K con soglia Bekker (0.60) e Functional Mask
  2. fig5b_wt_analysisMD_observables.png (Osservabili da analysisMD.tcl del Prof. Gautieri):
     - A: Contenuto di Alpha-Helix (%) nel tempo (300 K vs 400 K)
     - B: Numero di Legami a Idrogeno (Hbonds) nel tempo (300 K vs 400 K)
     - C: Raggio di Girazione Rgyr (Å) nel tempo (300 K vs 400 K)
     - D: SASA - Solvent Accessible Surface Area (Å^2) + Q-value per Repeat ANK1..4 a 400 K
"""

import os
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("06_presentation", exist_ok=True)

# Caricamento dati principali
rmsd_300  = np.loadtxt("02_wt_baseline/rmsd_300K.dat")
rmsd_400  = np.loadtxt("02_wt_baseline/rmsd_400K.dat")
rmsf_300  = np.loadtxt("02_wt_baseline/rmsf_300K.dat")
rmsf_400  = np.loadtxt("02_wt_baseline/rmsf_400K.dat")
q_300     = np.loadtxt("02_wt_baseline/q_time_300K.dat")
q_400     = np.loadtxt("02_wt_baseline/q_time_400K.dat")
sasa_300  = np.loadtxt("02_wt_baseline/sasa_300K.dat")
sasa_400  = np.loadtxt("02_wt_baseline/sasa_400K.dat")
rgyr_300  = np.loadtxt("02_wt_baseline/rgyr_300K.dat")
rgyr_400  = np.loadtxt("02_wt_baseline/rgyr_400K.dat")
alpha_300 = np.loadtxt("02_wt_baseline/alpha_300K.dat")
alpha_400 = np.loadtxt("02_wt_baseline/alpha_400K.dat")
hb_300    = np.loadtxt("02_wt_baseline/hbonds_300K.dat")
hb_400    = np.loadtxt("02_wt_baseline/hbonds_400K.dat")

# Caricamento q_per_residue_400K.dat
res_ids, n_cont, mean_qi, is_mask = [], [], [], []
res_names = []
with open("02_wt_baseline/q_per_residue_400K.dat", "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        res_ids.append(int(parts[0]))
        res_names.append(parts[1])
        n_cont.append(int(parts[2]))
        mean_qi.append(float(parts[3]))
        is_mask.append(int(parts[4]))

res_ids = np.array(res_ids)
n_cont  = np.array(n_cont)
mean_qi = np.array(mean_qi)
is_mask = np.array(is_mask)

def add_ank_bands(ax):
    bands = [
        (1, 34, "#e8f1fa", "ANK1"),
        (35, 67, "#e6f7f5", "ANK2"),
        (68, 100, "#fff2e6", "ANK3"),
        (101, 118, "#eafaea", "ANK4")
    ]
    for r_start, r_end, col, lbl in bands:
        ax.axvspan(r_start, r_end, color=col, alpha=0.65, zorder=0)
        ax.axvline(r_end, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)

# ==============================================================================
# FIGURA 5A: Pannello Principale (RMSD, Q-value, RMSF, <Q_i> per residuo)
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 9.5), dpi=300)
plt.subplots_adjust(hspace=0.34, wspace=0.24)

# 1. Pannello A: RMSD(t)
ax = axes[0, 0]
ax.plot(rmsd_300[:, 0], rmsd_300[:, 1], color="#1f77b4", lw=1.8, label="WT 300 K (Controllo)")
ax.plot(rmsd_400[:, 0], rmsd_400[:, 1], color="#d62728", lw=1.8, label="WT 400 K (Stress Termico)")
ax.set_title("A. Backbone RMSD vs Frame", fontsize=12, fontweight="bold")
ax.set_xlabel("Trajectory Frame")
ax.set_ylabel("Backbone RMSD (Å)")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(frameon=True)

# 2. Pannello B: Q-value nel tempo (Globale e per Repeat a 400 K)
ax = axes[0, 1]
ax.plot(q_300[:, 0], q_300[:, 1], color="#1f77b4", lw=2.0, label="300 K - Q: all-all")
ax.plot(q_400[:, 0], q_400[:, 1], color="#d62728", lw=2.2, label="400 K - Q: all-all")
ax.plot(q_400[:, 0], q_400[:, 2], color="#ff7f0e", lw=1.5, ls="--", label="400 K - Q: hydrophilic-all")
ax.set_title("B. Fraction of Native Contacts $Q(t)$", fontsize=12, fontweight="bold")
ax.set_xlabel("Trajectory Frame")
ax.set_ylabel("Native Contacts Fraction $Q$")
ax.set_ylim(0.4, 1.03)
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(frameon=True)

# 3. Pannello C: RMSF per residuo (300 K vs 400 K)
ax = axes[1, 0]
add_ank_bands(ax)
ax.plot(rmsf_300[:, 0], rmsf_300[:, 1], color="#1f77b4", lw=1.8, label="WT 300 K")
ax.plot(rmsf_400[:, 0], rmsf_400[:, 1], color="#d62728", lw=1.8, label="WT 400 K")
ax.set_title("C. Per-Residue Flexibility (RMSF $C_\\alpha$) across Ankyrin Repeats", fontsize=12, fontweight="bold")
ax.set_xlabel("Residue Number (1..118)")
ax.set_ylabel("RMSF (Å)")
ax.set_xlim(1, 118)
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(frameon=True)

# 4. Pannello D: <Q_i> per singolo residuo a 400 K
ax = axes[1, 1]
add_ank_bands(ax)
ax.plot(res_ids, mean_qi, color="#2b2b2b", lw=1.5, zorder=2, label="$\\langle Q_i \\rangle$ at 400 K")
cand_mask = (is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (mean_qi < 0.65)
func_pts  = (is_mask == 1)
ax.scatter(res_ids[func_pts], mean_qi[func_pts], color="#1f77b4", s=28, zorder=3, label="Functional Mask (CP interface)")
ax.scatter(res_ids[cand_mask], mean_qi[cand_mask], color="#d62728", s=48, zorder=4, label="Unstable Candidates (Scaffold)")
ax.axhline(0.60, color="#d62728", linestyle="--", lw=1.4, label="Bekker Threshold ($Q_i = 0.60$)")

non_masked_idx = np.where((is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 20))[0]
if len(non_masked_idx) > 0:
    worst_idx = non_masked_idx[np.argsort(mean_qi[non_masked_idx])[:5]]
    for idx in worst_idx:
        ax.annotate(
            f"{res_names[idx]}{res_ids[idx]}",
            (res_ids[idx], mean_qi[idx]),
            textcoords="offset points",
            xytext=(0, -13),
            ha="center",
            fontsize=8,
            fontweight="bold",
            color="#b30000"
        )

ax.set_title("D. Per-Residue Native Contact Retention $\\langle Q_i \\rangle$ at 400 K", fontsize=12, fontweight="bold")
ax.set_xlabel("Residue Number (1..118)")
ax.set_ylabel("Mean Contact Fraction $\\langle Q_i \\rangle$")
ax.set_xlim(1, 118)
ax.set_ylim(0.0, 1.05)
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(loc="lower left", fontsize=8.5, frameon=True)

out_pdf1 = "06_presentation/fig5_wt_baseline_plots.pdf"
plt.savefig(out_pdf1, format="pdf", bbox_inches="tight")
plt.close(fig)
print(f"[OK] Grafico principale vettoriale salvato in: {out_pdf1}")

# ==============================================================================
# FIGURA 5B: Osservabili da analysisMD.tcl (Alpha-Helix %, Hbonds, Rgyr, Q per Repeat)
# ==============================================================================
fig2, axes2 = plt.subplots(2, 2, figsize=(14, 9.5), dpi=300)
plt.subplots_adjust(hspace=0.34, wspace=0.24)

# A: Alpha-Helix (%)
ax = axes2[0, 0]
ax.plot(alpha_300[:, 0], alpha_300[:, 1] * 100.0, color="#1f77b4", lw=1.8, label="WT 300 K")
ax.plot(alpha_400[:, 0], alpha_400[:, 1] * 100.0, color="#d62728", lw=1.8, label="WT 400 K")
ax.set_title("A. $\\alpha$-Helix Content (%) vs Frame (STRIDE)", fontsize=12, fontweight="bold")
ax.set_xlabel("Trajectory Frame")
ax.set_ylabel("$\\alpha$-Helicity (%)")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(frameon=True)

# B: Hydrogen Bonds + Rgyr
ax = axes2[0, 1]
ax.plot(hb_300[:, 0], hb_300[:, 1], color="#1f77b4", lw=1.8, label="WT 300 K - H-bonds")
ax.plot(hb_400[:, 0], hb_400[:, 1], color="#d62728", lw=1.8, label="WT 400 K - H-bonds")
ax.set_title("B. Intra-Protein Hydrogen Bonds vs Frame", fontsize=12, fontweight="bold")
ax.set_xlabel("Trajectory Frame")
ax.set_ylabel("Number of H-Bonds")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(frameon=True)

# C: Rgyr & SASA
ax = axes2[1, 0]
ax.plot(rgyr_300[:, 0], rgyr_300[:, 1], color="#1f77b4", lw=1.8, label="WT 300 K - $R_g$")
ax.plot(rgyr_400[:, 0], rgyr_400[:, 1], color="#d62728", lw=1.8, label="WT 400 K - $R_g$")
ax.set_title("C. Radius of Gyration ($R_g$) vs Frame", fontsize=12, fontweight="bold")
ax.set_xlabel("Trajectory Frame")
ax.set_ylabel("$R_g$ (Å)")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(frameon=True)

# D: Q-value per singolo Ankyrin Repeat (ANK1..ANK4) a 400 K
ax = axes2[1, 1]
ax.plot(q_400[:, 0], q_400[:, 4], lw=1.8, color="#1f77b4", label="ANK1 (N-cap, res 1-34)")
ax.plot(q_400[:, 0], q_400[:, 5], lw=1.8, color="#17becf", label="ANK2 (Core, res 35-67)")
ax.plot(q_400[:, 0], q_400[:, 6], lw=1.8, color="#ff7f0e", label="ANK3 (Core, res 68-100)")
ax.plot(q_400[:, 0], q_400[:, 7], lw=1.8, color="#2ca02c", label="ANK4 (C-cap, res 101-118)")
ax.set_title("D. Repeat-Resolved $Q(t)$ at 400 K (ANK1–ANK4)", fontsize=12, fontweight="bold")
ax.set_xlabel("Trajectory Frame")
ax.set_ylabel("Native Contacts Fraction $Q$")
ax.set_ylim(0.4, 1.03)
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(frameon=True)

out_pdf2 = "06_presentation/fig5b_wt_analysisMD_observables.pdf"
plt.savefig(out_pdf2, format="pdf", bbox_inches="tight")
plt.close(fig2)
print(f"[OK] Grafico osservabili analysisMD.tcl vettoriale salvato in: {out_pdf2}")
