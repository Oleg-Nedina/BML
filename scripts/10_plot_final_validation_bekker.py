#!/usr/bin/env python3
"""
10_plot_final_validation_bekker.py
Punto 5 della Scaletta del Prof. Gautieri (Laboratory_Project_2026.pdf, Slide 3, 4, 11-13, 15, 22-23)
& Livelli 3, 4, 5 della Politica di Validazione (ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md):

1. Legge da 05_analysis/ tutte le repliche multi-seed analizzate da 09_analyze_matched_validation.tcl
   per WT e varianti mutate (K37L, T103L, K37L_T103L, H68G, K37L_H68G_T103L).
2. Calcola sulla finestra finale equilibrata (ultimi 30% dei frame, come in Bekker et al. 2019):
   - Media ± Deviazione Standard tra le repliche (μ ± σ) per:
     * Q: hydrophilic-all (r=0.84 con Tm!), Q: all-all (r=0.79), Q: hydrophobic-small (r=0.76)
     * Q risolto per Ankyrin Repeat (Q_ANK1..4) e per Interfaccia (Q_int12, Q_int23, Q_int34)
     * Numero di residui dello scaffold localmente instabili (<Q_i>_400K < 0.60, N_cont >= 30)
     * Backbone RMSD, Functional Mask 3AAA RMSD, Rgyr, SASA, % Alpha-Helix (STRIDE), H-Bonds
   - p-value del test t di Welch rispetto al Wild-Type e Ciclo Epistatico Dinamico del Doppio Mutante.
3. Genera:
   - 05_analysis/bekker_table2_and_epistasis_report.txt (Table II completa da slide!)
   - 06_presentation/fig9_matched_q_time_and_repeats.pdf
   - 06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf
   - 06_presentation/fig11_epistatic_cycle_and_observables.pdf
"""

import os
import glob
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

os.makedirs("05_analysis", exist_ok=True)
os.makedirs("06_presentation", exist_ok=True)

VARIANT_ORDER = ["WT", "K37L", "T103L", "K37L_T103L", "H68G", "K37L_H68G_T103L"]
VARIANT_COLORS = {
    "WT":              "#d62728",
    "K37L":            "#1f77b4",
    "T103L":           "#ff7f0e",
    "K37L_T103L":      "#2ca02c",
    "H68G":            "#17becf",
    "K37L_H68G_T103L": "#9467bd",
}

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

        res_ids, n_cont, qi_vals, is_mask = [], [], [], []
        if os.path.exists(qi_f):
            with open(qi_f) as f:
                for line in f:
                    if not line.startswith("#") and line.strip():
                        p = line.split()
                        res_ids.append(int(p[0]))
                        n_cont.append(int(p[2]))
                        qi_vals.append(float(p[3]))
                        is_mask.append(int(p[4]))
        res_ids = np.array(res_ids)
        n_cont = np.array(n_cont)
        qi_vals = np.array(qi_vals)
        is_mask = np.array(is_mask)

        # Conteggio residui dello scaffold interno localmente instabili (<Qi> < 0.60, N_0 >= 30, no terminali 1-3, 116-118)
        scaffold_unstable = int(np.sum((is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.60)))
        total_unstable = int(np.sum((res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.60)))

        nf = len(q_arr)
        f_start = max(1, int(nf * 0.70))
        q_win = np.mean(q_arr[f_start:, 1:], axis=0)  # Q_all, Q_phil, Q_ps, Q_ANK1..4, Q_int12..34
        obs_win = np.mean(obs_arr[f_start:, 1:], axis=0) if obs_arr is not None else np.zeros(6)
        hb_win = float(np.mean(hb_arr[f_start:, 1])) if hb_arr is not None else 0.0

        reps.append({
            "q_arr": q_arr,
            "obs_arr": obs_arr,
            "rmsf_arr": rmsf_arr,
            "qi_vals": qi_vals,
            "res_ids": res_ids,
            "is_mask": is_mask,
            "n_cont": n_cont,
            "q_win": q_win,
            "obs_win": obs_win,
            "hb_win": hb_win,
            "scaffold_unstable": scaffold_unstable,
            "total_unstable": total_unstable,
        })

    data_by_var[vname] = reps

if not data_by_var:
    print("ATTENZIONE: Nessun file trovato in 05_analysis/. Esegui prima: vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl")
    raise SystemExit(0)

# Scrittura della Table II di Bekker et al. (2019) e Report Statistico in 05_analysis/bekker_table2_and_epistasis_report.txt
report_out = "05_analysis/bekker_table2_and_epistasis_report.txt"
with open(report_out, "w") as f:
    f.write("========================================================================================================================================================\n")
    f.write(" TABLE II (BEKKER, MA & KAMIYA 2019, p. 434): STATISTICHE DI STABILITA' A 400 K RISPETTO AL WILD-TYPE (FINESTRA EQUILIBRATA FINALE 30%)\n")
    f.write("========================================================================================================================================================\n")
    f.write(f"{'Sistema':<17} {'Seeds':<6} {'RMSD_BB (Å)':<15} {'Q: phil-all (r=.84)':<21} {'Q: all-all (r=.79)':<20} {'Q: phob-small':<16} {'N_Unstable(Scaf)':<18} {'p-val(vs WT)'}\n")
    f.write("-" * 135 + "\n")

    wt_q_all_samples = [r["q_win"][0] for r in data_by_var.get("WT", [])]

    for vname, reps in data_by_var.items():
        n_s = len(reps)
        rmsd_vals = [r["obs_win"][0] for r in reps]
        q_all_vals = [r["q_win"][0] for r in reps]
        q_phil_vals = [r["q_win"][1] for r in reps]
        q_ps_vals = [r["q_win"][2] for r in reps]
        unst_vals = [r["scaffold_unstable"] for r in reps]

        if vname != "WT" and len(wt_q_all_samples) >= 2 and len(q_all_vals) >= 2:
            _, pval = stats.ttest_ind(q_all_vals, wt_q_all_samples, equal_var=False)
            p_str = f"p = {pval:.4f}"
        else:
            p_str = "   -   "

        f.write(
            f"{vname:<17s} {n_s:<6d} "
            f"{np.mean(rmsd_vals):.2f} ± {np.std(rmsd_vals):.2f}    "
            f"{np.mean(q_phil_vals):.4f} ± {np.std(q_phil_vals):.4f}      "
            f"{np.mean(q_all_vals):.4f} ± {np.std(q_all_vals):.4f}     "
            f"{np.mean(q_ps_vals):.4f} ± {np.std(q_ps_vals):.4f}   "
            f"{np.mean(unst_vals):.1f} ± {np.std(unst_vals):.1f}          "
            f"{p_str}\n"
        )

    f.write("\n========================================================================================================================================================\n")
    f.write(" STABILITA' RISOLTA PER MODULO ANKYRIN REPEAT (Q_ANK1..4), INTERFACCE INTER-REPEAT E OSSERVABILI analysisMD.tcl\n")
    f.write("========================================================================================================================================================\n")
    f.write(f"{'Sistema':<17} {'Q_ANK1':<12} {'Q_ANK2':<12} {'Q_ANK3':<12} {'Q_ANK4':<12} {'Q_int(1-2/2-3/3-4)':<24} {'Alpha-Helix(%)':<16} {'Rgyr (Å)':<14} {'Mask3AAA_RMSD'}\n")
    f.write("-" * 145 + "\n")
    for vname, reps in data_by_var.items():
        qa1 = np.mean([r["q_win"][3] for r in reps])
        qa2 = np.mean([r["q_win"][4] for r in reps])
        qa3 = np.mean([r["q_win"][5] for r in reps])
        qa4 = np.mean([r["q_win"][6] for r in reps])
        qi12 = np.mean([r["q_win"][7] for r in reps])
        qi23 = np.mean([r["q_win"][8] for r in reps])
        qi34 = np.mean([r["q_win"][9] for r in reps])
        alph = np.mean([r["obs_win"][4] * 100.0 for r in reps])
        rgyr = np.mean([r["obs_win"][2] for r in reps])
        m_rmsd = np.mean([r["obs_win"][1] for r in reps])
        f.write(
            f"{vname:<17s} {qa1:<12.4f} {qa2:<12.4f} {qa3:<12.4f} {qa4:<12.4f} "
            f"{qi12:.3f} / {qi23:.3f} / {qi34:<8.3f} {alph:<16.1f} {rgyr:<14.2f} {m_rmsd:.2f} Å\n"
        )

print(f"[OK] Report Bekker Table II salvato in: {report_out}")
