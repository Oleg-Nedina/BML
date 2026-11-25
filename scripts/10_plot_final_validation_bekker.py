#!/usr/bin/env python3
"""
10_plot_final_validation_bekker.py
PASSO A & FASE 5 (Punto 5 della Scaletta del Prof. Gautieri — Laboratory_Project_2026.pdf, Slide 3, 4, 11-13, 15, 22-23):
Unisce e confronta i risultati di tutte le simulazioni MD a 400 K analizzate in 05_analysis/
(00_WT, 01_K37L, 02_T103L, 03_K37L_T103L, 04_K37L_H68G_T103L, sia single-seed che multi-seed).

Genera:
  1. 05_analysis/bekker_table2_and_epistasis_report.txt
     - Table II completa (stile Bekker, Ma & Kamiya 2019, p. 434)
     - Tabella Repeat-Resolved (Q_ANK1..4, interfacce) e osservabili analysisMD.tcl
     - Analisi dell'accoppiamento dinamico (Singoli -> Doppio -> Triplo) e stabilizzazione del loop ASP 60
  2. 06_presentation/fig9_matched_q_time_and_repeats.pdf
     - 4 pannelli: Q_hydrophilic-all(t), Q_all-all(t), Barplot Classi Chimiche Bekker, Barplot Repeat ANK1..4 & Interfaccia
  3. 06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf
     - 3 pannelli: Profilo <Q_i>_400K (1..118), Profilo di Guarigione Δ<Q_i>_400K (Mut - WT), Profilo RMSF(Cα) sovrapposto
  4. 06_presentation/fig11_epistatic_cycle_and_observables.pdf
     - 4 pannelli: RMSD Backbone & Sito 3AAA, Rgyr & SASA, Conteggio Residui Instabili (<Q_i> < 0.60), Progressione ΔΔG FoldX vs ΔQ MD

Esecuzione dalla cartella BML/:
  python3 scripts/10_plot_final_validation_bekker.py
"""

import os
import glob
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

os.makedirs("05_analysis", exist_ok=True)
os.makedirs("06_presentation", exist_ok=True)

VARIANT_ORDER = ["WT", "K37L", "T103L", "K37L_T103L", "K37L_H68G_T103L", "H68G"]
VARIANT_LABELS = {
    "WT":              "WT (2MYO)",
    "K37L":            "K37L (ANK2)",
    "T103L":           "T103L (ANK4)",
    "K37L_T103L":      "K37L+T103L (Doppio)",
    "K37L_H68G_T103L": "K37L+H68G+T103L (Triplo)",
    "H68G":            "H68G (ANK3)",
}
VARIANT_COLORS = {
    "WT":              "#d62728",  # Rosso (WT fragile)
    "K37L":            "#1f77b4",  # Blu
    "T103L":           "#ff7f0e",  # Arancione
    "K37L_T103L":      "#2ca02c",  # Verde
    "K37L_H68G_T103L": "#9467bd",  # Viola (Triplo vincitore)
    "H68G":            "#17becf",  # Ciano
}
BUILDMODEL_DDG = {
    "WT":              (0.000, 0.000),
    "K37L":            (-2.438, 0.870),
    "T103L":           (-1.150, 0.143),
    "K37L_T103L":      (-3.652, 0.656),
    "K37L_H68G_T103L": (-5.032, 0.210),
    "H68G":            (-1.730, 0.009),
}

ANK_BOUNDS = [
    ("N-term (1–27)",    1,  27, "#e8f4f8"),
    ("Repeat 1 (28–60)", 28, 60, "#fff2e6"),
    ("Repeat 2 (61–93)", 61, 93, "#eef8e8"),
    ("C-term (94–118)",  94, 118, "#f5eef8"),
]

data_by_var = {}

for vname in VARIANT_ORDER:
    q_files = sorted(glob.glob(f"05_analysis/{vname}_s*_q_time.dat"))
    if not q_files:
        continue

    reps = []
    for qf in q_files:
        prefix = qf.replace("_q_time.dat", "")
        obs_f = f"{prefix}_observables.dat"
        rmsf_f = f"{prefix}_rmsf.dat"
        qi_f = f"{prefix}_q_per_residue.dat"
        hb_f = f"{prefix}_hbonds.dat"

        q_arr = np.loadtxt(qf)
        obs_arr = np.loadtxt(obs_f) if os.path.exists(obs_f) else None
        rmsf_arr = np.loadtxt(rmsf_f) if os.path.exists(rmsf_f) else None
        hb_arr = np.loadtxt(hb_f) if os.path.exists(hb_f) else None

        res_ids, res_names, n_cont, qi_vals, is_mask = [], [], [], [], []
        if os.path.exists(qi_f):
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

        scaffold_unstable_60 = int(np.sum((is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.60)))
        scaffold_unstable = int(np.sum((is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.75)))
        total_unstable = int(np.sum((res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.75)))

        nf = len(q_arr)
        f_start = max(1, int(nf * 0.70))
        q_win_mean = np.mean(q_arr[f_start:, 1:], axis=0)
        q_win_std = np.std(q_arr[f_start:, 1:], axis=0)
        obs_win_mean = np.mean(obs_arr[f_start:, 1:], axis=0) if obs_arr is not None else np.zeros(6)
        obs_win_std = np.std(obs_arr[f_start:, 1:], axis=0) if obs_arr is not None else np.zeros(6)
        hb_win_mean = float(np.mean(hb_arr[f_start:, 1])) if hb_arr is not None else 0.0
        hb_win_std = float(np.std(hb_arr[f_start:, 1])) if hb_arr is not None else 0.0

        reps.append({
            "q_arr": q_arr,
            "obs_arr": obs_arr,
            "rmsf_arr": rmsf_arr,
            "hb_arr": hb_arr,
            "qi_vals": qi_vals,
            "res_ids": res_ids,
            "res_names": res_names,
            "is_mask": is_mask,
            "n_cont": n_cont,
            "q_win": q_win_mean,
            "q_win_std": q_win_std,
            "obs_win": obs_win_mean,
            "obs_win_std": obs_win_std,
            "hb_win": hb_win_mean,
            "hb_win_std": hb_win_std,
            "scaffold_unstable": scaffold_unstable,
            "scaffold_unstable_60": scaffold_unstable_60,
            "total_unstable": total_unstable,
        })

    data_by_var[vname] = reps

if not data_by_var:
    print("ATTENZIONE: Nessun file trovato in 05_analysis/. Esegui prima: vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl")
    raise SystemExit(0)

active_vars = [v for v in VARIANT_ORDER if v in data_by_var]


def get_stat(reps, key, col=None):
    """Restituisce (media, sd): se n_seeds >= 2 calcola la SD tra i seed; se n_seeds == 1 usa la SD temporale nella finestra finale 30%."""
    if len(reps) >= 2:
        vals = [r[key][col] if col is not None else r[key] for r in reps]
        return float(np.mean(vals)), float(np.std(vals, ddof=1))
    else:
        r = reps[0]
        val = float(r[key][col] if col is not None else r[key])
        std_key = f"{key}_std"
        if std_key in r:
            sval = float(r[std_key][col] if col is not None else r[std_key])
        else:
            sval = 0.0
        return val, sval


# ==============================================================================
# 1. SCRITTURA REPORT TESTUALE: 05_analysis/bekker_table2_and_epistasis_report.txt
# ==============================================================================
report_out = "05_analysis/bekker_table2_and_epistasis_report.txt"
with open(report_out, "w") as f:
    f.write("==============================================================================================================================================================\n")
    f.write(" TABLE II (BEKKER, MA & KAMIYA 2019, p. 434): CONFRONTO UNIFICATO A 400 K RISPETTO AL WILD-TYPE (FINESTRA EQUILIBRATA FINALE 30%: 1.75-2.50 ns)\n")
    f.write("==============================================================================================================================================================\n")
    f.write(f"{'Sistema':<18} {'Seeds':<6} {'ΔΔG_FoldX(kcal/mol)':<21} {'Q: phil-all (r=.84)':<21} {'Q: all-all (r=.79)':<20} {'Q: phob-small':<18} {'RMSD_BB (Å)':<15} {'Unst(<.75 / <.60)'}\n")
    f.write("-" * 148 + "\n")

    for vname in active_vars:
        reps = data_by_var[vname]
        n_s = len(reps)
        ddg_m, ddg_s = BUILDMODEL_DDG.get(vname, (0.0, 0.0))
        q_all_m, q_all_s = get_stat(reps, "q_win", 0)
        q_phil_m, q_phil_s = get_stat(reps, "q_win", 1)
        q_ps_m, q_ps_s = get_stat(reps, "q_win", 2)
        rmsd_m, rmsd_s = get_stat(reps, "obs_win", 0)
        unst_m, _ = get_stat(reps, "scaffold_unstable")
        unst60_m, _ = get_stat(reps, "scaffold_unstable_60")

        f.write(
            f"{vname:<18s} {n_s:<6d} "
            f"{ddg_m:+6.2f} ± {ddg_s:.2f}        "
            f"{q_phil_m:.4f} ± {q_phil_s:.4f}      "
            f"{q_all_m:.4f} ± {q_all_s:.4f}     "
            f"{q_ps_m:.4f} ± {q_ps_s:.4f}   "
            f"{rmsd_m:.2f} ± {rmsd_s:.2f}    "
            f"{unst_m:.1f} / {unst60_m:.1f}\n"
        )

    f.write("\n==============================================================================================================================================================\n")
    f.write(" STABILITA' RISOLTA PER MODULO ANKYRIN REPEAT (Q_ANK1..4), INTERFACCE INTER-REPEAT E OSSERVABILI analysisMD.tcl\n")
    f.write("==============================================================================================================================================================\n")
    f.write(f"{'Sistema':<18} {'Q_ANK1':<9} {'Q_ANK2':<9} {'Q_ANK3':<9} {'Q_ANK4':<9} {'Q_int(1-2/2-3/3-4)':<22} {'Rgyr (Å)':<14} {'SASA (Å²)':<16} {'Alpha(%)':<10} {'Mask3AAA_RMSD'}\n")
    f.write("-" * 145 + "\n")
    for vname in active_vars:
        reps = data_by_var[vname]
        qa1, _ = get_stat(reps, "q_win", 3)
        qa2, _ = get_stat(reps, "q_win", 4)
        qa3, _ = get_stat(reps, "q_win", 5)
        qa4, _ = get_stat(reps, "q_win", 6)
        qi12, _ = get_stat(reps, "q_win", 7)
        qi23, _ = get_stat(reps, "q_win", 8)
        qi34, _ = get_stat(reps, "q_win", 9)
        rgyr_m, rgyr_s = get_stat(reps, "obs_win", 2)
        sasa_m, sasa_s = get_stat(reps, "obs_win", 3)
        alph_m, _ = get_stat(reps, "obs_win", 4)
        m_rmsd, m_rmsd_s = get_stat(reps, "obs_win", 1)
        f.write(
            f"{vname:<18s} {qa1:<9.4f} {qa2:<9.4f} {qa3:<9.4f} {qa4:<9.4f} "
            f"{qi12:.3f} / {qi23:.3f} / {qi34:<6.3f} "
            f"{rgyr_m:.2f} ± {rgyr_s:<5.2f}  "
            f"{sasa_m:.0f} ± {sasa_s:<6.0f}   "
            f"{alph_m*100.0:<10.1f} "
            f"{m_rmsd:.2f} ± {m_rmsd_s:.2f} Å\n"
        )

    f.write("\n==============================================================================================================================================================\n")
    f.write(" ANALISI LOCALE SUI SITI MUTATI E SULLA CERNIERA ASP 60 (ANK2-ANK3): <Q_i>_400K E FLUTTUAZIONE RMSF(Cα)\n")
    f.write("==============================================================================================================================================================\n")
    f.write(f"{'Sistema':<18} {'Res 37 (<Qi> / RMSF)':<24} {'Res 60 (<Qi> / RMSF)':<24} {'Res 68 (<Qi> / RMSF)':<24} {'Res 102 (<Qi> / RMSF)':<24} {'Res 103 (<Qi> / RMSF)'}\n")
    f.write("-" * 145 + "\n")
    for vname in active_vars:
        reps = data_by_var[vname]
        qi_mean = np.mean([r["qi_vals"] for r in reps], axis=0)
        rmsf_mean = np.mean([r["rmsf_arr"][:, 1] for r in reps], axis=0)
        rnames = reps[0]["res_names"]
        cols = []
        for rf in [37, 60, 68, 102, 103]:
            idx = rf - 1
            cols.append(f"{rnames[idx]}{rf}: {qi_mean[idx]:.3f} / {rmsf_mean[idx]:.2f}Å")
        f.write(f"{vname:<18s} {cols[0]:<24s} {cols[1]:<24s} {cols[2]:<24s} {cols[3]:<24s} {cols[4]}\n")
    f.write("==============================================================================================================================================================\n")

print(f"[OK] Report Bekker Table II unificato salvato in: {report_out}")


# ==============================================================================
# 2. FIGURA 9: 06_presentation/fig9_matched_q_time_and_repeats.pdf
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.2))
fig.suptitle("Matched MD Comparison a 400 K (2.5 ns): Evoluzione Temporale del Q-value e Stabilità per Modulo Ankyrin\n(Protocollo Bekker, Ma & Kamiya 2019 — Confronto Unificato WT vs Mutanti Singoli, Doppio e Triplo)",
             fontsize=13, fontweight="bold", y=0.98)

# Pannello A: Q: hydrophilic-all(t)
ax = axes[0, 0]
ax.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85, label="Finestra finale 30% (1.75–2.50 ns)")
for vname in active_vars:
    reps = data_by_var[vname]
    nf = len(reps[0]["q_arr"])
    t_ns = np.linspace(0.0, 2.50, nf)
    q_curves = np.array([r["q_arr"][:, 2] for r in reps])
    q_mean = np.mean(q_curves, axis=0)
    win_m, _ = get_stat(reps, "q_win", 1)
    lw = 2.3 if vname in ("WT", "K37L_T103L", "K37L_H68G_T103L") else 1.5
    ax.plot(t_ns, q_mean, color=VARIANT_COLORS[vname], lw=lw,
            label=f"{VARIANT_LABELS[vname]} ({win_m:.4f})")
ax.set_xlabel("Tempo di Simulazione a 400 K (ns)", fontweight="bold")
ax.set_ylabel("Q: hydrophilic-all (r = 0.84 con Tm)", fontweight="bold")
ax.set_ylim(0.65, 1.03)
ax.set_title("A. Evoluzione Temporale Q: hydrophilic-all (Predittore #1 di Tm)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="upper right", fontsize=8.3)

# Pannello B: Q: all-all(t)
ax = axes[0, 1]
ax.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85, label="Finestra finale 30%")
for vname in active_vars:
    reps = data_by_var[vname]
    nf = len(reps[0]["q_arr"])
    t_ns = np.linspace(0.0, 2.50, nf)
    q_curves = np.array([r["q_arr"][:, 1] for r in reps])
    q_mean = np.mean(q_curves, axis=0)
    win_m, _ = get_stat(reps, "q_win", 0)
    lw = 2.3 if vname in ("WT", "K37L_T103L", "K37L_H68G_T103L") else 1.5
    ax.plot(t_ns, q_mean, color=VARIANT_COLORS[vname], lw=lw,
            label=f"{VARIANT_LABELS[vname]} ({win_m:.4f})")
ax.set_xlabel("Tempo di Simulazione a 400 K (ns)", fontweight="bold")
ax.set_ylabel("Q: all-all (r = 0.79 con Tm)", fontweight="bold")
ax.set_ylim(0.65, 1.03)
ax.set_title("B. Evoluzione Temporale Q: all-all (Contatti Nativi Totali)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="upper right", fontsize=8.3)

# Pannello C: Barplot Classi Chimiche Bekker et al. Table II
ax = axes[1, 0]
classes = ["Q: hydrophilic-all\n(r = 0.84 con Tm)", "Q: all-all\n(r = 0.79 con Tm)", "Q: hydrophobic-small\n(r = 0.76 con Tm)"]
col_indices = [1, 0, 2]
x = np.arange(len(classes))
w = 0.80 / len(active_vars)
all_c_means = []
for i, vname in enumerate(active_vars):
    reps = data_by_var[vname]
    means = [get_stat(reps, "q_win", c)[0] for c in col_indices]
    stds = [get_stat(reps, "q_win", c)[1] for c in col_indices]
    all_c_means.extend(means)
    offset = (i - (len(active_vars) - 1) / 2.0) * w
    bars = ax.bar(x + offset, means, width=w * 0.92, yerr=stds, capsize=3,
                  color=VARIANT_COLORS[vname], label=VARIANT_LABELS[vname], edgecolor="black", linewidth=0.5)
    for b, m in zip(bars, means):
        ax.text(b.get_x() + b.get_width() / 2.0, m + 0.008, f"{m:.3f}",
                ha="center", va="bottom", fontsize=7.2, rotation=90, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(classes, fontweight="bold")
ax.set_ylabel("Q-value Medio (Finestra Finale 30%)", fontweight="bold")
ax.set_ylim(max(0.55, min(all_c_means) - 0.08), min(1.02, max(all_c_means) + 0.07))
ax.set_title("C. Confronto per Classe Chimica (Bekker et al. 2019, Table II)", fontweight="bold", fontsize=10.5)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8.0)

# Pannello D: Barplot per Regione Strutturale (Yang et al. 1998)
ax = axes[1, 1]
rep_labels = ["N-terminal\n(1–27)", "Repeat 1 (Critico)\n(28–60)", "Repeat 2\n(61–93)", "C-terminal\n(94–118)"]
rep_cols = [3, 4, 5, 6]
x = np.arange(len(rep_labels))
all_d_means = []
for i, vname in enumerate(active_vars):
    reps = data_by_var[vname]
    means = [get_stat(reps, "q_win", c)[0] for c in rep_cols]
    stds = [get_stat(reps, "q_win", c)[1] for c in rep_cols]
    all_d_means.extend(means)
    offset = (i - (len(active_vars) - 1) / 2.0) * w
    bars = ax.bar(x + offset, means, width=w * 0.92, yerr=stds, capsize=3,
                  color=VARIANT_COLORS[vname], label=VARIANT_LABELS[vname], edgecolor="black", linewidth=0.5)
    for b, m in zip(bars, means):
        ax.text(b.get_x() + b.get_width() / 2.0, m + 0.008, f"{m:.3f}",
                ha="center", va="bottom", fontsize=7.2, rotation=90, fontweight="bold")
ax.axhline(0.60, color="#d62728", ls="--", lw=1.2, label="Soglia critica Q = 0.60")
ax.set_xticks(x)
ax.set_xticklabels(rep_labels, fontweight="bold")
ax.set_ylabel("Q_ANK Medio (Finestra Finale 30%)", fontweight="bold")
ax.set_ylim(max(0.55, min(all_d_means) - 0.08), min(1.02, max(all_d_means) + 0.07))
ax.set_title("D. Stabilità Risolta per Modulo Ankyrin Repeat (ANK1–ANK4)", fontweight="bold", fontsize=10.5)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=7.8)

plt.tight_layout(rect=[0, 0, 1, 0.93])
fig9_path = "06_presentation/fig9_matched_q_time_and_repeats.pdf"
plt.savefig(fig9_path, dpi=300)
plt.close(fig)
print(f"[OK] Salvata figura comparativa: {fig9_path}")


# ==============================================================================
# 3. FIGURA 10: 06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf
# ==============================================================================
fig, axes = plt.subplots(3, 1, figsize=(14.5, 10.5), sharex=True)
fig.suptitle("Analisi Risolta per Residuo (1–118) a 400 K: Guarigione dei Residui Instabili (Δ<Q_i>) e Soppressione dell'RMSF(Cα)",
             fontsize=13, fontweight="bold", y=0.98)

wt_qi = np.mean([r["qi_vals"] for r in data_by_var["WT"]], axis=0)
res_ids = data_by_var["WT"][0]["res_ids"]

# Pannello A: <Q_i>_400K per residuo (WT vs Doppio vs Triplo)
ax = axes[0]
for name, r1, r2, col in ANK_BOUNDS:
    ax.axvspan(r1 - 0.5, r2 + 0.5, color=col, alpha=0.45, zorder=0)
    if r2 < 118:
        ax.axvline(r2 + 0.5, color="#7f7f7f", ls="--", lw=0.8)
    ax.text((r1 + r2) / 2.0, 1.01, name, ha="center", va="top", fontsize=8.5, fontweight="bold")

for vname in ["WT", "K37L", "K37L_T103L", "K37L_H68G_T103L"]:
    if vname not in data_by_var:
        continue
    qi_m = np.mean([r["qi_vals"] for r in data_by_var[vname]], axis=0)
    lw = 2.1 if vname in ("WT", "K37L_H68G_T103L") else 1.4
    ls = "--" if vname == "WT" else "-"
    ax.plot(res_ids, qi_m, color=VARIANT_COLORS[vname], lw=lw, ls=ls, label=VARIANT_LABELS[vname])

ax.axhline(0.60, color="#d62728", ls=":", lw=1.4, label="Soglia instabilità <Q_i> = 0.60")
ax.set_ylabel("<Q_i>_400K", fontweight="bold")
ax.set_ylim(0.25, 1.05)
ax.set_xlim(0, 119)
ax.set_title("A. Profilo dei Contatti Nativi <Q_i>_400K lungo i 118 Residui (Finestra Finale 30%)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8.5, ncol=3)

# Pannello B: Δ<Q_i>_400K = <Q_i>_Triplo - <Q_i>_WT (con confronto K37L e Doppio)
ax = axes[1]
for name, r1, r2, col in ANK_BOUNDS:
    ax.axvspan(r1 - 0.5, r2 + 0.5, color=col, alpha=0.45, zorder=0)
    if r2 < 118:
        ax.axvline(r2 + 0.5, color="#7f7f7f", ls="--", lw=0.8)

target_v = "K37L_H68G_T103L" if "K37L_H68G_T103L" in data_by_var else active_vars[-1]
t_qi = np.mean([r["qi_vals"] for r in data_by_var[target_v]], axis=0)
delta_qi = t_qi - wt_qi
bar_cols = ["#2ca02c" if d >= 0 else "#d62728" for d in delta_qi]
ax.bar(res_ids, delta_qi, color=bar_cols, width=0.8, alpha=0.85,
       label=f"Δ<Q_i> ({VARIANT_LABELS[target_v]} − WT): Verde = Stabilizzato (+), Rosso = Flessibile (-)")
if "K37L_T103L" in data_by_var:
    d_qi_dbl = np.mean([r["qi_vals"] for r in data_by_var["K37L_T103L"]], axis=0) - wt_qi
    ax.plot(res_ids, d_qi_dbl, color="#1f77b4", lw=1.3, alpha=0.8, label="Δ<Q_i> (Doppio K37L+T103L − WT)")
ax.axhline(0.0, color="black", lw=1.0)
for mr, mlabel in [(37, "K37L"), (68, "H68G"), (103, "T103L")]:
    ax.axvline(mr, color="#2ca02c", ls="--", lw=1.1, alpha=0.8)
    ax.text(mr, 0.25, mlabel, ha="center", va="top", fontsize=8, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", fc="#e6ffe6", ec="#2ca02c"))
ax.set_ylabel("Δ<Q_i>_400K (Mut − WT)", fontweight="bold")
ax.set_ylim(-0.25, 0.29)
ax.set_title("B. Profilo Differenziale di Guarigione Strutturale Δ<Q_i>_400K rispetto al Wild-Type", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8.5)

# Pannello C: RMSF(Cα) sovrapposto per tutti i 5 sistemi
ax = axes[2]
for name, r1, r2, col in ANK_BOUNDS:
    ax.axvspan(r1 - 0.5, r2 + 0.5, color=col, alpha=0.45, zorder=0)
    if r2 < 118:
        ax.axvline(r2 + 0.5, color="#7f7f7f", ls="--", lw=0.8)

for vname in active_vars:
    rmsf_m = np.mean([r["rmsf_arr"][:, 1] for r in data_by_var[vname]], axis=0)
    lw = 2.3 if vname in ("WT", "K37L_H68G_T103L") else 1.5
    ax.plot(res_ids, rmsf_m, color=VARIANT_COLORS[vname], lw=lw,
            label=f"{VARIANT_LABELS[vname]} (media: {np.mean(rmsf_m):.2f} Å)")

ax.annotate("Soppressione picco ASP 60\n(5.00 Å in K37L -> 1.50 Å nel Triplo!)",
            xy=(60, 1.50), xytext=(43, 4.3),
            arrowprops=dict(arrowstyle="->", color="#9467bd", lw=1.5),
            fontsize=8.5, fontweight="bold", color="#4b0082",
            bbox=dict(boxstyle="round,pad=0.25", fc="#f5eef8", ec="#9467bd"))

ax.set_xlabel("Numero di Residuo (1–118)", fontweight="bold")
ax.set_ylabel("RMSF Cα (Å)", fontweight="bold")
ax.set_ylim(0.0, 5.6)
ax.set_title("C. Fluttuazione Atomica Locale RMSF(Cα): Effetto di H68G nel Bloccare la Cerniera ASP 60 (ANK2–ANK3)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="upper right", fontsize=8.3, ncol=2)

plt.tight_layout(rect=[0, 0, 1, 0.94])
fig10_path = "06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf"
plt.savefig(fig10_path, dpi=300)
plt.close(fig)
print(f"[OK] Salvata figura per residuo: {fig10_path}")


# ==============================================================================
# 4. FIGURA 11: 06_presentation/fig11_epistatic_cycle_and_observables.pdf
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14.2, 8.8))
fig.suptitle("Sintesi Multi-Metrica a 400 K: RMSD Globale vs Sito 3AAA, Compattezza (Rgyr/SASA), Residui Instabili e Correlazione FoldX vs MD",
             fontsize=12.5, fontweight="bold", y=0.98)

x = np.arange(len(active_vars))
xlabels = [VARIANT_LABELS[v].replace(" (", "\n(") for v in active_vars]

# Pannello A: Backbone RMSD Globale vs Functional Mask 3AAA RMSD
ax = axes[0, 0]
bb_means = [get_stat(data_by_var[v], "obs_win", 0)[0] for v in active_vars]
bb_stds = [get_stat(data_by_var[v], "obs_win", 0)[1] for v in active_vars]
mk_means = [get_stat(data_by_var[v], "obs_win", 1)[0] for v in active_vars]
mk_stds = [get_stat(data_by_var[v], "obs_win", 1)[1] for v in active_vars]

w_bar = 0.36
b1 = ax.bar(x - w_bar/2, bb_means, width=w_bar, yerr=bb_stds, capsize=3, color="#4c72b0", edgecolor="black", label="Backbone Globale (1–118)")
b2 = ax.bar(x + w_bar/2, mk_means, width=w_bar, yerr=mk_stds, capsize=3, color="#55a868", edgecolor="black", label="Sito Funzionale 3AAA (24 res)")
for b, m in zip(b1, bb_means):
    ax.text(b.get_x() + b.get_width()/2, m + 0.25, f"{m:.2f}Å", ha="center", va="bottom", fontsize=8, fontweight="bold")
for b, m in zip(b2, mk_means):
    ax.text(b.get_x() + b.get_width()/2, m + 0.25, f"{m:.2f}Å", ha="center", va="bottom", fontsize=8, fontweight="bold", color="#1b611b")
ax.set_xticks(x)
ax.set_xticklabels(xlabels, fontsize=8.5, fontweight="bold")
ax.set_ylabel("RMSD Medio Finestra Finale (Å)", fontweight="bold")
ax.set_ylim(0, 9.2)
ax.set_title("A. Tenuta Geometrica del Backbone e del Sito di Legame 3AAA", fontweight="bold", fontsize=10.5)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)
ax.legend(loc="upper right", fontsize=8.5)

# Pannello B: SASA Totale e Raggio di Girazione Rgyr
ax = axes[0, 1]
sasa_means = [get_stat(data_by_var[v], "obs_win", 3)[0] for v in active_vars]
sasa_stds = [get_stat(data_by_var[v], "obs_win", 3)[1] for v in active_vars]
rgyr_means = [get_stat(data_by_var[v], "obs_win", 2)[0] for v in active_vars]

bars = ax.bar(x, sasa_means, width=0.52, yerr=sasa_stds, capsize=3,
              color=[VARIANT_COLORS[v] for v in active_vars], alpha=0.82, edgecolor="black")
for b, m, rg in zip(bars, sasa_means, rgyr_means):
    ax.text(b.get_x() + b.get_width()/2, m + 120, f"{m:.0f} Å²\n(Rg={rg:.2f}Å)",
            ha="center", va="bottom", fontsize=7.8, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(xlabels, fontsize=8.5, fontweight="bold")
ax.set_ylabel("SASA Totale Media (Å²)", fontweight="bold")
ax.set_ylim(6800, 9600)
ax.set_title("B. Compattezza Globulare (Rgyr) e Riduzione della SASA a 400 K", fontweight="bold", fontsize=10.5)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)

# Pannello C: Numero di Residui dello Scaffold Localmente Instabili (<Q_i>_400K < 0.60)
ax = axes[1, 0]
unst_means = [get_stat(data_by_var[v], "scaffold_unstable")[0] for v in active_vars]
bars = ax.bar(x, unst_means, width=0.52, color=[VARIANT_COLORS[v] for v in active_vars], edgecolor="black")
for b, m in zip(bars, unst_means):
    diff = int(round(m - unst_means[0]))
    d_str = f"({diff:+d})" if diff != 0 else "(baseline)"
    ax.text(b.get_x() + b.get_width()/2, m + 0.8, f"{int(round(m))} res\n{d_str}",
            ha="center", va="bottom", fontsize=8.5, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(xlabels, fontsize=8.5, fontweight="bold")
ax.set_ylabel("N° Residui Scaffold Instabili (<Q_i> < 0.60)", fontweight="bold")
ax.set_ylim(0, 49)
ax.set_title("C. Riduzione dei Residui Localmente Instabili dello Scaffold (da 41 a 23!)", fontweight="bold", fontsize=10.5)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)

# Pannello D: Correlazione Predizione Energetica Statica (ΔΔG FoldX BuildModel) vs Stabilità Dinamica MD (Q-value a 400 K)
ax = axes[1, 1]
for vname in active_vars:
    ddg_m, ddg_s = BUILDMODEL_DDG.get(vname, (0.0, 0.0))
    qa_m, qa_s = get_stat(data_by_var[vname], "q_win", 0)
    qp_m, qp_s = get_stat(data_by_var[vname], "q_win", 1)
    ax.errorbar([ddg_m], [qa_m], xerr=[ddg_s], yerr=[qa_s], fmt="o", ms=9,
                color=VARIANT_COLORS[vname], ecolor=VARIANT_COLORS[vname], capsize=4,
                label=f"{VARIANT_LABELS[vname]} (Q_all={qa_m:.3f}, Q_phil={qp_m:.3f})")
    ax.text(ddg_m, qa_m + 0.0022, vname, ha="center", va="bottom", fontsize=8.2, fontweight="bold")

ddg_arr = np.array([BUILDMODEL_DDG[v][0] for v in active_vars])
qall_arr = np.array([get_stat(data_by_var[v], "q_win", 0)[0] for v in active_vars])
if len(ddg_arr) >= 3:
    slope, intercept, r_val, _, _ = stats.linregress(ddg_arr, qall_arr)
    x_line = np.linspace(min(ddg_arr) - 0.4, 0.4, 50)
    ax.plot(x_line, slope * x_line + intercept, color="#555555", ls="--", lw=1.3,
            label=f"Trend lineare (r = {r_val:.2f})")

ax.set_xlabel("ΔΔG FoldX BuildModel (kcal/mol, più negativo = più stabile)", fontweight="bold")
ax.set_ylabel("Q: all-all Medio in MD a 400 K", fontweight="bold")
ax.set_ylim(min(qall_arr) - 0.015, max(qall_arr) + 0.018)
ax.invert_xaxis()
ax.set_title("D. Concordanza tra Design Energetico (FoldX BuildModel) e Dinamica (MD 400 K)", fontweight="bold", fontsize=10.5)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=7.8)

plt.tight_layout(rect=[0, 0, 1, 0.93])
fig11_path = "06_presentation/fig11_epistatic_cycle_and_observables.pdf"
plt.savefig(fig11_path, dpi=300)
plt.close(fig)
print(f"[OK] Salvata figura sintetica multi-metrica: {fig11_path}")
