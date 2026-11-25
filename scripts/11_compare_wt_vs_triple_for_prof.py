#!/usr/bin/env python3
"""
11_compare_wt_vs_triple_for_prof.py
Genera una cartella dedicata pronta per il colloquio con il Prof. Gautieri:
  06_presentation/presentazione_prof_WT_vs_Triplo/

Contenuto generato (confronto diretto 1-vs-1: Wild-Type vs Triplo Mutante K37L+H68G+T103L a 400 K):
  1. 00_TABELLA_CONFRONTO_WT_vs_TRIPLO.txt      (Tabella numerica completa WT vs Triplo + Δ%)
  2. 01_WT_vs_Triplo_Q_value_tempo_e_moduli.pdf (+ .png)
  3. 02_WT_vs_Triplo_RMSD_Rgyr_SASA_3AAA.pdf    (+ .png)
  4. 03_WT_vs_Triplo_per_residuo_Qi_e_RMSF.pdf  (+ .png)

Esecuzione dalla cartella BML/:
  python3 scripts/11_compare_wt_vs_triple_for_prof.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

OUT_DIR = "06_presentation/presentazione_prof_WT_vs_Triplo"
os.makedirs(OUT_DIR, exist_ok=True)

WT_PREFIX = "05_analysis/WT_s1"
TR_PREFIX = "05_analysis/K37L_H68G_T103L_s1"

COL_WT = "#d62728"  # Rosso
COL_TR = "#6a3d9a"  # Viola scuro (Triplo vincitore)

ANK_BOUNDS = [
    ("N-term (1–27)",    1,  27, "#e8f4f8"),
    ("Repeat 1 (28–60)", 28, 60, "#fff2e6"),
    ("Repeat 2 (61–93)", 61, 93, "#eef8e8"),
    ("C-term (94–118)",  94, 118, "#f5eef8"),
]


def load_system(prefix):
    qf = f"{prefix}_q_time.dat"
    obs_f = f"{prefix}_observables.dat"
    rmsf_f = f"{prefix}_rmsf.dat"
    qi_f = f"{prefix}_q_per_residue.dat"
    hb_f = f"{prefix}_hbonds.dat"

    for fp in (qf, obs_f, rmsf_f, qi_f):
        if not os.path.exists(fp):
            raise FileNotFoundError(
                f"File mancante: {fp}. Esegui prima: vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl"
            )

    q_arr = np.loadtxt(qf)
    obs_arr = np.loadtxt(obs_f)
    rmsf_arr = np.loadtxt(rmsf_f)
    hb_arr = np.loadtxt(hb_f) if os.path.exists(hb_f) else None

    res_ids, res_names, n_cont, qi_vals, is_mask = [], [], [], [], []
    with open(qi_f) as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                p = line.split()
                res_ids.append(int(p[0]))
                res_names.append(p[1])
                n_cont.append(int(p[2]))
                qi_vals.append(float(p[3]))
                is_mask.append(int(p[4]))

    res_ids = np.array(res_ids)
    n_cont = np.array(n_cont)
    qi_vals = np.array(qi_vals)
    is_mask = np.array(is_mask)

    nf = len(q_arr)
    f_start = max(1, int(nf * 0.70))
    t_ns = np.linspace(0.0, 2.50, nf)

    unst_75 = int(np.sum((is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.75)))
    unst_60 = int(np.sum((is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.60)))

    return {
        "q": q_arr,
        "obs": obs_arr,
        "rmsf": rmsf_arr,
        "hb": hb_arr,
        "res_ids": res_ids,
        "res_names": res_names,
        "n_cont": n_cont,
        "qi": qi_vals,
        "is_mask": is_mask,
        "nf": nf,
        "f_start": f_start,
        "t_ns": t_ns,
        "q_mean": np.mean(q_arr[f_start:, 1:], axis=0),
        "q_std": np.std(q_arr[f_start:, 1:], axis=0),
        "obs_mean": np.mean(obs_arr[f_start:, 1:], axis=0),
        "obs_std": np.std(obs_arr[f_start:, 1:], axis=0),
        "hb_mean": float(np.mean(hb_arr[f_start:, 1])) if hb_arr is not None else 0.0,
        "hb_std": float(np.std(hb_arr[f_start:, 1])) if hb_arr is not None else 0.0,
        "unst_75": unst_75,
        "unst_60": unst_60,
    }


wt = load_system(WT_PREFIX)
tr = load_system(TR_PREFIX)


def save_both(fig, stem):
    pdf_p = os.path.join(OUT_DIR, f"{stem}.pdf")
    png_p = os.path.join(OUT_DIR, f"{stem}.png")
    fig.savefig(pdf_p, dpi=300, bbox_inches="tight")
    fig.savefig(png_p, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Salvato: {pdf_p} (+ .png)")


# ==============================================================================
# 0. TABELLA COMPARATIVA TESTUALE 1-VS-1: 00_TABELLA_CONFRONTO_WT_vs_TRIPLO.txt
# ==============================================================================
tab_path = os.path.join(OUT_DIR, "00_TABELLA_CONFRONTO_WT_vs_TRIPLO.txt")
with open(tab_path, "w") as f:
    f.write("=====================================================================================================\n")
    f.write(" CONFRONTO DIRETTO PER IL PROF. GAUTIERI: WILD-TYPE (2MYO) vs TRIPLO MUTANTE (K37L + H68G + T103L)\n")
    f.write(" Condizioni Matched: CHARMM36m + TIP3P (12 Å) + 0.15 M NaCl | 400 K | 2.50 ns | Finestra Finale 30%\n")
    f.write(" Riferimento Q-value: Frame 0 della traiettoria DCD minimizzata -> Q(0) = 1.0000\n")
    f.write("=====================================================================================================\n")
    f.write(f"{'Osservabile (1.75–2.50 ns a 400 K)':<42} {'Wild-Type (WT)':<22} {'Triplo (K37L+H68G+T103L)':<25} {'Delta (Triplo - WT)'}\n")
    f.write("-" * 101 + "\n")

    rows = [
        ("ΔΔG FoldX BuildModel (n=5 run, kcal/mol)", "+0.00 ± 0.00", "-5.03 ± 0.21", "-5.03 kcal/mol"),
        ("Q: all-all (Bekker r=0.79 con Tm)",
         f"{wt['q_mean'][0]:.4f} ± {wt['q_std'][0]:.4f}",
         f"{tr['q_mean'][0]:.4f} ± {tr['q_std'][0]:.4f}",
         f"{tr['q_mean'][0] - wt['q_mean'][0]:+.4f} ({(tr['q_mean'][0]-wt['q_mean'][0])*100:+.2f}%)"),
        ("Q: hydrophilic-all (Bekker r=0.84 con Tm)",
         f"{wt['q_mean'][1]:.4f} ± {wt['q_std'][1]:.4f}",
         f"{tr['q_mean'][1]:.4f} ± {tr['q_std'][1]:.4f}",
         f"{tr['q_mean'][1] - wt['q_mean'][1]:+.4f} ({(tr['q_mean'][1]-wt['q_mean'][1])*100:+.2f}%)"),
        ("Q: hydrophobic-small (Bekker r=0.76)",
         f"{wt['q_mean'][2]:.4f} ± {wt['q_std'][2]:.4f}",
         f"{tr['q_mean'][2]:.4f} ± {tr['q_std'][2]:.4f}",
         f"{tr['q_mean'][2] - wt['q_mean'][2]:+.4f} ({(tr['q_mean'][2]-wt['q_mean'][2])*100:+.2f}%)"),
        ("Q_Nterm (N-terminal, res 1–27)", f"{wt['q_mean'][3]:.4f}", f"{tr['q_mean'][3]:.4f}", f"{tr['q_mean'][3]-wt['q_mean'][3]:+.4f}"),
        ("Q_Rep1 (Repeat 1 critico, res 28–60)", f"{wt['q_mean'][4]:.4f}", f"{tr['q_mean'][4]:.4f}", f"{tr['q_mean'][4]-wt['q_mean'][4]:+.4f}"),
        ("Q_Rep2 (Repeat 2, res 61–93)", f"{wt['q_mean'][5]:.4f}", f"{tr['q_mean'][5]:.4f}", f"{tr['q_mean'][5]-wt['q_mean'][5]:+.4f}"),
        ("Q_Cterm (C-terminal, res 94–118)", f"{wt['q_mean'][6]:.4f}", f"{tr['q_mean'][6]:.4f}", f"{tr['q_mean'][6]-wt['q_mean'][6]:+.4f}"),
        ("Interfaccia Q_int(N-term – Rep1)", f"{wt['q_mean'][7]:.4f}", f"{tr['q_mean'][7]:.4f}", f"{tr['q_mean'][7]-wt['q_mean'][7]:+.4f}"),
        ("Backbone RMSD Globale (Å)",
         f"{wt['obs_mean'][0]:.2f} ± {wt['obs_std'][0]:.2f} Å",
         f"{tr['obs_mean'][0]:.2f} ± {tr['obs_std'][0]:.2f} Å",
         f"{tr['obs_mean'][0]-wt['obs_mean'][0]:+.2f} Å ({(tr['obs_mean'][0]-wt['obs_mean'][0])/wt['obs_mean'][0]*100:+.1f}%)"),
        ("Sito Funzionale 3AAA RMSD (24 res, Å)",
         f"{wt['obs_mean'][1]:.2f} ± {wt['obs_std'][1]:.2f} Å",
         f"{tr['obs_mean'][1]:.2f} ± {tr['obs_std'][1]:.2f} Å",
         f"{tr['obs_mean'][1]-wt['obs_mean'][1]:+.2f} Å ({(tr['obs_mean'][1]-wt['obs_mean'][1])/wt['obs_mean'][1]*100:+.1f}%)"),
        ("Raggio di Girazione Rgyr (Å)",
         f"{wt['obs_mean'][2]:.2f} ± {wt['obs_std'][2]:.2f} Å",
         f"{tr['obs_mean'][2]:.2f} ± {tr['obs_std'][2]:.2f} Å",
         f"{tr['obs_mean'][2]-wt['obs_mean'][2]:+.2f} Å"),
        ("SASA Totale (Å²)",
         f"{wt['obs_mean'][3]:.1f} ± {wt['obs_std'][3]:.1f} Å²",
         f"{tr['obs_mean'][3]:.1f} ± {tr['obs_std'][3]:.1f} Å²",
         f"{tr['obs_mean'][3]-wt['obs_mean'][3]:+.1f} Å²"),
        ("% Alpha-Helix (STRIDE)", f"{wt['obs_mean'][4]*100:.2f}%", f"{tr['obs_mean'][4]*100:.2f}%", f"{(tr['obs_mean'][4]-wt['obs_mean'][4])*100:+.2f}%"),
        ("H-Bonds Intra-Proteici Unici", f"{wt['hb_mean']:.1f} ± {wt['hb_std']:.1f}", f"{tr['hb_mean']:.1f} ± {tr['hb_std']:.1f}", f"{tr['hb_mean']-wt['hb_mean']:+.1f}"),
        ("Residui Scaffold Deboli (<Q_i> < 0.75)", f"{wt['unst_75']} residui", f"{tr['unst_75']} residui", f"{tr['unst_75']-wt['unst_75']:+d} residui"),
        ("Residui Scaffold Critici (<Q_i> < 0.60)", f"{wt['unst_60']} residui", f"{tr['unst_60']} residui", f"{tr['unst_60']-wt['unst_60']:+d} residui"),
    ]
    for name, v_wt, v_tr, diff in rows:
        f.write(f"{name:<42s} {v_wt:<22s} {v_tr:<25s} {diff}\n")

    f.write("-----------------------------------------------------------------------------------------------------\n")
    f.write(" FLUTTUAZIONE LOCALE RMSF(Cα) E <Q_i>_400K SUI SITI CHIAVE:\n")
    for rf in [37, 60, 68, 102, 103]:
        idx = rf - 1
        f.write(
            f"  - Res {rf:3d} ({wt['res_names'][idx]} -> {tr['res_names'][idx]}): "
            f"RMSF {wt['rmsf'][idx,1]:.2f} Å -> {tr['rmsf'][idx,1]:.2f} Å ({tr['rmsf'][idx,1]-wt['rmsf'][idx,1]:+.2f} Å) | "
            f"<Q_i> {wt['qi'][idx]:.3f} -> {tr['qi'][idx]:.3f} ({tr['qi'][idx]-wt['qi'][idx]:+.3f})\n"
        )
    f.write("=====================================================================================================\n")

print(f"[OK] Tabella comparativa salvata in: {tab_path}")


# ==============================================================================
# 1. PLOT 1: 01_WT_vs_Triplo_Q_value_tempo_e_moduli
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14.0, 8.8))
fig.suptitle("Wild-Type (2MYO) vs Triplo Mutante (K37L + H68G + T103L) a 400 K\nFrazione di Contatti Nativi Q(t) (da Q(0) = 1.00) per Classi Chimiche e Moduli Ankyrin",
             fontsize=13, fontweight="bold", y=0.98)

# A: Q: all-all(t)
ax = axes[0, 0]
ax.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85, label="Finestra finale 30% (1.75–2.50 ns)")
ax.plot(wt["t_ns"], wt["q"][:, 1], color=COL_WT, lw=1.9, label=f"Wild-Type (media fin: {wt['q_mean'][0]:.4f})")
ax.plot(tr["t_ns"], tr["q"][:, 1], color=COL_TR, lw=2.1, label=f"Triplo K37L+H68G+T103L (media fin: {tr['q_mean'][0]:.4f})")
ax.set_xlabel("Tempo di Simulazione a 400 K (ns)", fontweight="bold")
ax.set_ylabel("Q: all-all", fontweight="bold")
ax.set_ylim(0.65, 1.03)
ax.set_title("A. Contatti Nativi Globali Q: all-all(t)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower left", fontsize=8.8)

# B: Q: hydrophilic-all(t)
ax = axes[0, 1]
ax.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85, label="Finestra finale 30%")
ax.plot(wt["t_ns"], wt["q"][:, 2], color=COL_WT, lw=1.9, label=f"Wild-Type (media fin: {wt['q_mean'][1]:.4f})")
ax.plot(tr["t_ns"], tr["q"][:, 2], color=COL_TR, lw=2.1, label=f"Triplo K37L+H68G+T103L (media fin: {tr['q_mean'][1]:.4f})")
ax.set_xlabel("Tempo di Simulazione a 400 K (ns)", fontweight="bold")
ax.set_ylabel("Q: hydrophilic-all (r = 0.84 con Tm)", fontweight="bold")
ax.set_ylim(0.65, 1.03)
ax.set_title("B. Contatti Polari Q: hydrophilic-all(t) (Bekker et al.)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower left", fontsize=8.8)

# C: Barplot Classi Chimiche Bekker
ax = axes[1, 0]
lbls = ["Q: hydrophilic-all\n(r=0.84 con Tm)", "Q: all-all\n(r=0.79 con Tm)", "Q: hydrophobic-small\n(r=0.76 con Tm)"]
idxs = [1, 0, 2]
x = np.arange(len(lbls))
w = 0.34
m_wt = [wt["q_mean"][i] for i in idxs]
s_wt = [wt["q_std"][i] for i in idxs]
m_tr = [tr["q_mean"][i] for i in idxs]
s_tr = [tr["q_std"][i] for i in idxs]
b1 = ax.bar(x - w/2, m_wt, width=w, yerr=s_wt, capsize=4, color=COL_WT, edgecolor="black", label="Wild-Type (WT)")
b2 = ax.bar(x + w/2, m_tr, width=w, yerr=s_tr, capsize=4, color=COL_TR, edgecolor="black", label="Triplo (K37L+H68G+T103L)")
for b, m in zip(b1, m_wt):
    ax.text(b.get_x()+b.get_width()/2, m+0.012, f"{m:.3f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=COL_WT)
for b, m in zip(b2, m_tr):
    ax.text(b.get_x()+b.get_width()/2, m+0.012, f"{m:.3f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=COL_TR)
ax.set_xticks(x)
ax.set_xticklabels(lbls, fontweight="bold")
ax.set_ylabel("Q-value Medio (Finestra Finale 30%)", fontweight="bold")
ax.set_ylim(max(0.55, min(m_wt+m_tr)-0.08), min(1.02, max(m_wt+m_tr)+0.08))
ax.set_title("C. Confronto per Classe Chimica (Bekker Table II)", fontweight="bold", fontsize=10.5)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8.8)

# D: Barplot Regioni Strutturali (Yang et al. 1998)
ax = axes[1, 1]
rlbls = ["N-terminal\n(1–27)", "Repeat 1 (Critico)\n(28–60)", "Repeat 2\n(61–93)", "C-terminal\n(94–118)"]
ridxs = [3, 4, 5, 6]
x = np.arange(len(rlbls))
m_wt = [wt["q_mean"][i] for i in ridxs]
m_tr = [tr["q_mean"][i] for i in ridxs]
b1 = ax.bar(x - w/2, m_wt, width=w, color=COL_WT, edgecolor="black", label="Wild-Type (WT)")
b2 = ax.bar(x + w/2, m_tr, width=w, color=COL_TR, edgecolor="black", label="Triplo (K37L+H68G+T103L)")
for b, m in zip(b1, m_wt):
    ax.text(b.get_x()+b.get_width()/2, m+0.01, f"{m:.3f}", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color=COL_WT)
for b, m in zip(b2, m_tr):
    ax.text(b.get_x()+b.get_width()/2, m+0.01, f"{m:.3f}", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color=COL_TR)
ax.set_xticks(x)
ax.set_xticklabels(rlbls, fontweight="bold")
ax.set_ylabel("Q_ANK Medio (Finestra Finale 30%)", fontweight="bold")
ax.set_ylim(max(0.55, min(m_wt+m_tr)-0.08), min(1.02, max(m_wt+m_tr)+0.08))
ax.set_title("D. Stabilità Risolta per Modulo Ankyrin Repeat", fontweight="bold", fontsize=10.5)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8.8)

plt.tight_layout(rect=[0, 0, 1, 0.93])
save_both(fig, "01_WT_vs_Triplo_Q_value_tempo_e_moduli")


# ==============================================================================
# 2. PLOT 2: 02_WT_vs_Triplo_RMSD_Rgyr_SASA_3AAA
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14.0, 8.6), sharex=True)
fig.suptitle("Wild-Type (2MYO) vs Triplo Mutante (K37L + H68G + T103L) a 400 K\nOsservabili Strutturali Globali e Protezione del Sito Funzionale 3AAA",
             fontsize=13, fontweight="bold", y=0.98)

panels = [
    (axes[0, 0], 1, "Backbone RMSD (Å)", "A. Deviazione Backbone Globale (1–118)"),
    (axes[0, 1], 2, "Functional Mask 3AAA RMSD (Å)", "B. Tenuta Sito Funzionale 3AAA (24 residui CP-binding)"),
    (axes[1, 0], 3, "Raggio di Girazione Rgyr (Å)", "C. Compattezza Globulare (Rgyr)"),
    (axes[1, 1], 4, "SASA Totale (Å²)", "D. Esposizione al Solvente (SASA)"),
]

for ax, col_i, ylab, title in panels:
    ax.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85, label="Finestra finale 30%")
    m_w = wt["obs_mean"][col_i - 1]
    m_t = tr["obs_mean"][col_i - 1]
    fmt = "{:.0f}" if col_i == 4 else "{:.2f}"
    ax.plot(wt["t_ns"], wt["obs"][:, col_i], color=COL_WT, lw=1.8, label=f"Wild-Type (fin: {fmt.format(m_w)})")
    ax.plot(tr["t_ns"], tr["obs"][:, col_i], color=COL_TR, lw=2.0, label=f"Triplo Mutante (fin: {fmt.format(m_t)})")
    ax.set_ylabel(ylab, fontweight="bold")
    ax.set_title(title, fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper left", fontsize=8.8)

axes[1, 0].set_xlabel("Tempo di Simulazione a 400 K (ns)", fontweight="bold")
axes[1, 1].set_xlabel("Tempo di Simulazione a 400 K (ns)", fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.93])
save_both(fig, "02_WT_vs_Triplo_RMSD_Rgyr_SASA_3AAA")


# ==============================================================================
# 3. PLOT 3: 03_WT_vs_Triplo_per_residuo_Qi_e_RMSF
# ==============================================================================
fig, axes = plt.subplots(3, 1, figsize=(14.2, 10.2), sharex=True)
fig.suptitle("Analisi Locale Residuo-per-Residuo (1–118) a 400 K: Wild-Type vs Triplo Mutante (K37L + H68G + T103L)",
             fontsize=13, fontweight="bold", y=0.98)

res_ids = wt["res_ids"]

# A: <Q_i>_400K sovrapposto
ax = axes[0]
for name, r1, r2, col in ANK_BOUNDS:
    ax.axvspan(r1-0.5, r2+0.5, color=col, alpha=0.45, zorder=0)
    if r2 < 118:
        ax.axvline(r2+0.5, color="#7f7f7f", ls="--", lw=0.8)
    ax.text((r1+r2)/2.0, 1.02, name, ha="center", va="top", fontsize=8.5, fontweight="bold")
ax.plot(res_ids, wt["qi"], color=COL_WT, lw=1.7, ls="--", label=f"Wild-Type (<Q_i><0.75: {wt['unst_75']} res)")
ax.plot(res_ids, tr["qi"], color=COL_TR, lw=2.0, label=f"Triplo Mutante (<Q_i><0.75: {tr['unst_75']} res)")
ax.axhline(0.75, color="#ff7f0e", ls=":", lw=1.3, label="Soglia debolezza <Q_i> = 0.75")
ax.set_ylabel("<Q_i>_400K", fontweight="bold")
ax.set_ylim(0.30, 1.06)
ax.set_xlim(0, 119)
ax.set_title("A. Ritenzione Locale dei Contatti Nativi <Q_i>_400K (Finestra Finale 30%)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8.8)

# B: Δ<Q_i> = Triplo - WT
ax = axes[1]
for name, r1, r2, col in ANK_BOUNDS:
    ax.axvspan(r1-0.5, r2+0.5, color=col, alpha=0.45, zorder=0)
    if r2 < 118:
        ax.axvline(r2+0.5, color="#7f7f7f", ls="--", lw=0.8)
d_qi = tr["qi"] - wt["qi"]
bcols = ["#2ca02c" if v >= 0 else "#d62728" for v in d_qi]
ax.bar(res_ids, d_qi, color=bcols, width=0.82, alpha=0.85,
       label="Δ<Q_i> = Triplo − WT (Verde = Stabilizzato, Rosso = Ridotto)")
ax.axhline(0.0, color="black", lw=1.0)
for mr, mtxt in [(37, "K37L"), (68, "H68G"), (103, "T103L")]:
    ax.axvline(mr, color="#2ca02c", ls="--", lw=1.2)
    ax.text(mr, max(0.22, float(np.max(d_qi))*0.92), mtxt, ha="center", va="top", fontsize=8.5, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", fc="#e6ffe6", ec="#2ca02c"))
ax.set_ylabel("Δ<Q_i>_400K", fontweight="bold")
ax.set_title("B. Guadagno di Contatti Nativi per Residuo (Triplo − Wild-Type)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8.8)

# C: RMSF(Cα) sovrapposto
ax = axes[2]
for name, r1, r2, col in ANK_BOUNDS:
    ax.axvspan(r1-0.5, r2+0.5, color=col, alpha=0.45, zorder=0)
    if r2 < 118:
        ax.axvline(r2+0.5, color="#7f7f7f", ls="--", lw=0.8)
ax.plot(res_ids, wt["rmsf"][:, 1], color=COL_WT, lw=1.9, ls="--",
        label=f"Wild-Type (RMSF medio: {np.mean(wt['rmsf'][:,1]):.2f} Å)")
ax.plot(res_ids, tr["rmsf"][:, 1], color=COL_TR, lw=2.2,
        label=f"Triplo Mutante (RMSF medio: {np.mean(tr['rmsf'][:,1]):.2f} Å)")
for mr in [37, 60, 68, 103]:
    idx = mr - 1
    ax.scatter([mr], [tr["rmsf"][idx, 1]], color="#2ca02c", s=55, zorder=5)
ax.set_xlabel("Numero di Residuo (1–118)", fontweight="bold")
ax.set_ylabel("RMSF Cα (Å)", fontweight="bold")
ax.set_ylim(0.0, max(5.2, float(np.max(wt["rmsf"][:, 1])) * 1.12))
ax.set_title("C. Abbattimento della Fluttuazione Atomica Locale RMSF(Cα) in tutta la proteina", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="upper left", fontsize=8.8)

plt.tight_layout(rect=[0, 0, 1, 0.94])
save_both(fig, "03_WT_vs_Triplo_per_residuo_Qi_e_RMSF")

print(f"\n>>> PACCHETTO PRESENTAZIONE PROF PRONTO IN: {OUT_DIR}/")
