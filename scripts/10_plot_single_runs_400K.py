#!/usr/bin/env python3
"""
10_plot_single_runs_400K.py
Genera i plot e il report numerico SEPARATI per ogni singola simulazione a 400 K:
  - WT               (Wild-Type 2MYO, Seed 1)       -> 05_analysis/single_runs/00_WT/
  - K37L             (Simulazione 1: Singolo ANK2)  -> 05_analysis/single_runs/01_K37L/
  - T103L            (Simulazione 2: Singolo ANK4)  -> 05_analysis/single_runs/02_T103L/
  - K37L_T103L       (Simulazione 3: Doppio ANK2+4) -> 05_analysis/single_runs/03_K37L_T103L/
  - K37L_H68G_T103L  (Simulazione 4: Triplo ANK234) -> 05_analysis/single_runs/04_K37L_H68G_T103L/
  - (Opzionale: H68G se presente)

Per OGNI simulazione disponibile in 05_analysis/<VARIANT>_s1_*, crea nella sua cartella dedicata:
  1. <VARIANT>_01_dynamics_overview.pdf        (RMSD BB + Sito 3AAA, Rgyr + SASA, % Alpha-Helix, H-Bonds)
  2. <VARIANT>_02_qvalue_repeats_and_classes.pdf (Q per classe chimica Bekker + Q per modulo ANK1..4 e interfacce)
  3. <VARIANT>_03_per_residue_qi_and_rmsf.pdf  (<Q_i>_400K per residuo 1..118 + RMSF CA per residuo)
  4. <VARIANT>_summary_report.txt              (Report numerico dettagliato della singola run + lista residui instabili)

Esecuzione dalla cartella BML/:
  python3 scripts/10_plot_single_runs_400K.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

RUN_SPECS = [
    ("WT",              "00_WT",              "Wild-Type (2MYO) — 400 K (Seed 1)",                         [],             "#d62728"),
    ("K37L",            "01_K37L",            "Simulazione 1: Mutante Singolo K37L (Repeat 1) — 400 K",    [37],           "#1f77b4"),
    ("T103L",           "02_T103L",           "Simulazione 2: Mutante Singolo T103L (C-terminal) — 400 K", [103],          "#ff7f0e"),
    ("K37L_T103L",      "03_K37L_T103L",      "Simulazione 3: Doppio Mutante K37L+T103L — 400 K",          [37, 103],      "#2ca02c"),
    ("K37L_H68G_T103L", "04_K37L_H68G_T103L", "Simulazione 4: Triplo Mutante K37L+H68G+T103L — 400 K",     [37, 68, 103],  "#9467bd"),
    ("H68G",            "05_H68G",            "Simulazione Extra: Mutante Singolo H68G (Repeat 2) — 400 K",[68],           "#17becf"),
]

ANK_BOUNDS = [
    ("N-term (1–27)",    1,  27, "#e8f4f8"),
    ("Repeat 1 (28–60)", 28, 60, "#fff2e6"),
    ("Repeat 2 (61–93)", 61, 93, "#eef8e8"),
    ("C-term (94–118)",  94, 118, "#f5eef8"),
]


def add_ank_shading(ax, ymin, ymax, label_top=True):
    for name, r1, r2, col in ANK_BOUNDS:
        ax.axvspan(r1 - 0.5, r2 + 0.5, color=col, alpha=0.45, zorder=0)
        if r2 < 118:
            ax.axvline(r2 + 0.5, color="#7f7f7f", linestyle="--", linewidth=0.8, alpha=0.7)
        if label_top:
            ax.text((r1 + r2) / 2.0, ymax - 0.04 * (ymax - ymin), name,
                    ha="center", va="top", fontsize=8.5, fontweight="bold", color="#333333")


def process_single_run(vname, folder_name, title_label, mut_resids, main_color):
    prefix = f"05_analysis/{vname}_s1"
    qf = f"{prefix}_q_time.dat"
    obs_f = f"{prefix}_observables.dat"
    rmsf_f = f"{prefix}_rmsf.dat"
    qi_f = f"{prefix}_q_per_residue.dat"
    hb_f = f"{prefix}_hbonds.dat"

    if not (os.path.exists(qf) and os.path.exists(obs_f) and os.path.exists(rmsf_f) and os.path.exists(qi_f)):
        return False

    out_dir = os.path.join("05_analysis", "single_runs", folder_name)
    os.makedirs(out_dir, exist_ok=True)

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
    # Ogni frame di produzione (dcdfreq=2000 a dt=2fs) corrisponde a 4 ps = 0.004 ns (625 frame = 2.50 ns)
    time_ns = np.linspace(0.0, 2.50, nf)
    f_start = max(1, int(nf * 0.70))
    t_win_start = time_ns[f_start]

    # ==========================================================================
    # PLOT 1: <VARIANT>_01_dynamics_overview.pdf (4 pannelli temporali)
    # ==========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), sharex=True)
    fig.suptitle(f"{title_label}\nEvoluzione Strutturale Globale (2.5 ns, CHARMM36m + TIP3P)",
                 fontsize=13, fontweight="bold", y=0.98)

    # Pannello A: Backbone RMSD & 3AAA Functional Mask RMSD
    ax = axes[0, 0]
    ax.axvspan(t_win_start, 2.50, color="#fff8dc", alpha=0.8, label="Finestra finale 30% (1.75–2.50 ns)")
    ax.plot(time_ns, obs_arr[:, 1], color=main_color, lw=1.8, label=f"Backbone Globale (media fin: {np.mean(obs_arr[f_start:,1]):.2f} Å)")
    ax.plot(time_ns, obs_arr[:, 2], color="#2b5c8f", lw=1.5, ls="--", label=f"Sito Funzionale 3AAA (media fin: {np.mean(obs_arr[f_start:,2]):.2f} Å)")
    ax.set_ylabel("RMSD rispetto a Frame 0 (Å)", fontweight="bold")
    ax.set_title("A. Deviazione del Backbone (Globale vs Sito 3AAA)", fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper left", fontsize=8.5)

    # Pannello B: Raggio di Girazione (Rgyr) e SASA
    ax = axes[0, 1]
    ax.axvspan(t_win_start, 2.50, color="#fff8dc", alpha=0.8)
    ax.plot(time_ns, obs_arr[:, 3], color="#d95f02", lw=1.8, label=f"Rgyr (media fin: {np.mean(obs_arr[f_start:,3]):.2f} Å)")
    ax.set_ylabel("Raggio di Girazione Rgyr (Å)", fontweight="bold", color="#d95f02")
    ax.set_title("B. Compattezza Globulare (Rgyr) ed Esposizione al Solvente (SASA)", fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax2 = ax.twinx()
    ax2.plot(time_ns, obs_arr[:, 4], color="#7570b3", lw=1.3, alpha=0.75, label=f"SASA (media fin: {np.mean(obs_arr[f_start:,4]):.0f} Å²)")
    ax2.set_ylabel("SASA Totale (Å²)", fontweight="bold", color="#7570b3")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=8.5)

    # Pannello C: % Alpha-Helix (STRIDE)
    ax = axes[1, 0]
    ax.axvspan(t_win_start, 2.50, color="#fff8dc", alpha=0.8)
    alpha_perc = obs_arr[:, 5] * 100.0
    ax.plot(time_ns, alpha_perc, color="#1b9e77", lw=1.7, label=f"% α-Helix STRIDE (media fin: {np.mean(alpha_perc[f_start:]):.1f}%)")
    ax.axhline(45.8, color="#555555", ls=":", lw=1.1, label="Riferimento WT 300 K (45.8%)")
    ax.set_xlabel("Tempo di Simulazione (ns)", fontweight="bold")
    ax.set_ylabel("Contenuto Elicoidale (%)", fontweight="bold")
    ax.set_ylim(15, 65)
    ax.set_title("C. Tenuta della Struttura Secondaria α-Elicoidale (STRIDE)", fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="lower left", fontsize=8.5)

    # Pannello D: Legami a Idrogeno Intra-Proteici (H-Bonds)
    ax = axes[1, 1]
    ax.axvspan(t_win_start, 2.50, color="#fff8dc", alpha=0.8)
    if hb_arr is not None and len(hb_arr) == nf:
        ax.plot(time_ns, hb_arr[:, 1], color="#e7298a", lw=1.4, alpha=0.85,
                label=f"H-Bonds Intra-Proteici (media fin: {np.mean(hb_arr[f_start:,1]):.1f})")
    ax.set_xlabel("Tempo di Simulazione (ns)", fontweight="bold")
    ax.set_ylabel("Numero di H-Bonds Unici", fontweight="bold")
    ax.set_title("D. Rete di Legami a Idrogeno Intra-Proteici (-type unique)", fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper right", fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    f1_path = os.path.join(out_dir, f"{vname}_01_dynamics_overview.pdf")
    plt.savefig(f1_path, dpi=300)
    plt.close(fig)

    # ==========================================================================
    # PLOT 2: <VARIANT>_02_qvalue_repeats_and_classes.pdf (2 pannelli Q-value)
    # ==========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.4), sharex=True)
    fig.suptitle(f"{title_label}\nFrazione dei Contatti Nativi Q(t) per Classe Chimica (Bekker et al. 2019) e Modulo Ankyrin",
                 fontsize=12.5, fontweight="bold", y=0.99)

    # Pannello A: Classi Chimiche di Bekker et al. (2019) Table II
    ax = axes[0]
    ax.axvspan(t_win_start, 2.50, color="#fff8dc", alpha=0.8, label="Finestra 30% finale")
    ax.plot(time_ns, q_arr[:, 1], color="#222222", lw=2.0, label=f"Q: all-all (r=0.79) -> {np.mean(q_arr[f_start:,1]):.4f}")
    ax.plot(time_ns, q_arr[:, 2], color="#1f77b4", lw=1.7, label=f"Q: hydrophilic-all (r=0.84) -> {np.mean(q_arr[f_start:,2]):.4f}")
    ax.plot(time_ns, q_arr[:, 3], color="#ff7f0e", lw=1.7, label=f"Q: hydrophobic-small (r=0.76) -> {np.mean(q_arr[f_start:,3]):.4f}")
    ax.set_xlabel("Tempo di Simulazione (ns)", fontweight="bold")
    ax.set_ylabel("Frazione di Contatti Nativi Q(t)", fontweight="bold")
    ax.set_ylim(0.45, 1.02)
    ax.set_title("A. Q-value per Classe Chimica (Bekker et al. 2019, Table II)", fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="lower left", fontsize=8.5)

    # Pannello B: Q-value per Regione Strutturale (Yang et al. 1998)
    ax = axes[1]
    ax.axvspan(t_win_start, 2.50, color="#fff8dc", alpha=0.8)
    ax.plot(time_ns, q_arr[:, 4], color="#1f77b4", lw=1.7, label=f"N-terminal (1–27) -> {np.mean(q_arr[f_start:,4]):.4f}")
    ax.plot(time_ns, q_arr[:, 5], color="#d62728", lw=2.1, label=f"Ankyrin Repeat 1 (28–60) -> {np.mean(q_arr[f_start:,5]):.4f}")
    ax.plot(time_ns, q_arr[:, 6], color="#2ca02c", lw=1.7, label=f"Ankyrin Repeat 2 (61–93) -> {np.mean(q_arr[f_start:,6]):.4f}")
    ax.plot(time_ns, q_arr[:, 7], color="#9467bd", lw=1.9, label=f"C-terminal (94–118) -> {np.mean(q_arr[f_start:,7]):.4f}")
    ax.axhline(0.60, color="#b22222", ls="--", lw=1.1, label="Soglia critica Q = 0.60")
    ax.set_xlabel("Tempo di Simulazione (ns)", fontweight="bold")
    ax.set_ylabel("Q-value per Regione (Yang et al. 1998)", fontweight="bold")
    ax.set_ylim(0.45, 1.02)
    ax.set_title("B. Stabilità Risolta per Regione Strutturale (Yang et al. 1998)", fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="lower left", fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    f2_path = os.path.join(out_dir, f"{vname}_02_qvalue_repeats_and_classes.pdf")
    plt.savefig(f2_path, dpi=300)
    plt.close(fig)

    # ==========================================================================
    # PLOT 3: <VARIANT>_03_per_residue_qi_and_rmsf.pdf (<Q_i>_400K e RMSF)
    # ==========================================================================
    fig, axes = plt.subplots(2, 1, figsize=(14.0, 8.2), sharex=True)
    fig.suptitle(f"{title_label}\nProfilo Locale per Residuo: Contatti Nativi <Q_i>_400K e Fluttuazione RMSF(Cα)",
                 fontsize=12.5, fontweight="bold", y=0.98)

    # Identificazione residui instabili dello scaffold (<Q_i> < 0.75 e < 0.60, N_cont >= 30, non terminali 1-3, 116-118)
    unstable_scaf_mask = (is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.75)
    severe_scaf_mask = (is_mask == 0) & (res_ids > 3) & (res_ids < 116) & (n_cont >= 30) & (qi_vals < 0.60)
    n_unst_scaf = int(np.sum(unstable_scaf_mask))
    n_sev_scaf = int(np.sum(severe_scaf_mask))

    # Pannello A: <Q_i>_400K per residuo
    ax = axes[0]
    add_ank_shading(ax, 0.15, 1.05, label_top=True)
    bar_colors = []
    for r, m, u, sev in zip(res_ids, is_mask, unstable_scaf_mask, severe_scaf_mask):
        if r in mut_resids:
            bar_colors.append("#2ca02c")  # Verde brillante per i siti mutati
        elif m == 1:
            bar_colors.append("#aec7e8")  # Azzurro per la Functional Mask 3AAA
        elif sev:
            bar_colors.append("#b22222")  # Rosso scuro per residui fortemente instabili (<0.60)
        elif u:
            bar_colors.append("#ff7f0e")  # Arancio/Rosso per residui deboli (<0.75)
        else:
            bar_colors.append("#4c72b0")  # Blu standard per scaffold stabile (>=0.75)

    ax.bar(res_ids, qi_vals, color=bar_colors, width=0.82, edgecolor="none", zorder=3)
    ax.axhline(0.75, color="#ff7f0e", ls="--", lw=1.3, zorder=4,
               label=f"Soglia debolezza <Q_i> < 0.75 (Residui scaffold < 0.75: {n_unst_scaf})")
    ax.axhline(0.60, color="#b22222", ls=":", lw=1.3, zorder=4,
               label=f"Soglia critica <Q_i> < 0.60 (Residui scaffold < 0.60: {n_sev_scaf})")

    # Evidenzia i residui mutati
    for mr in mut_resids:
        idx = int(mr - 1)
        ax.scatter([mr], [qi_vals[idx] + 0.03], color="#2ca02c", s=70, marker="v", zorder=5)
        ax.text(mr, min(0.98, qi_vals[idx] + 0.06), f"Mut {res_names[idx]}{mr}\n({qi_vals[idx]:.2f})",
                ha="center", va="bottom", fontsize=8, fontweight="bold", color="#1b611b",
                bbox=dict(boxstyle="round,pad=0.2", fc="#e6ffe6", ec="#2ca02c", lw=0.8))

    ax.set_ylabel("<Q_i>_400K (ultimi 30% frame)", fontweight="bold")
    ax.set_ylim(0.15, 1.05)
    ax.set_xlim(0, 119)
    ax.set_title(f"A. Tenuta Locale dei Contatti Nativi <Q_i>_400K per Residuo (Rosso = Scaffold < 0.60 [{n_unst_scaf} residui]; Azzurro = Sito 3AAA)",
                 fontweight="bold", fontsize=10.5)
    ax.grid(True, axis="y", linestyle=":", alpha=0.5)
    ax.legend(loc="lower right", fontsize=8.5)

    # Pannello B: RMSF(CA) per residuo
    ax = axes[1]
    ymax_rmsf = max(4.5, float(np.max(rmsf_arr[:, 1])) * 1.12)
    add_ank_shading(ax, 0.0, ymax_rmsf, label_top=False)
    ax.plot(rmsf_arr[:, 0], rmsf_arr[:, 1], color=main_color, lw=2.0, marker="o", ms=3.0, zorder=3,
            label=f"RMSF Cα (media globale: {np.mean(rmsf_arr[:,1]):.2f} Å)")
    for mr in mut_resids:
        idx = int(mr - 1)
        ax.scatter([mr], [rmsf_arr[idx, 1]], color="#2ca02c", s=75, zorder=5,
                   label=f"Sito mutato {res_names[idx]}{mr} (RMSF={rmsf_arr[idx,1]:.2f} Å)")
    ax.set_xlabel("Numero di Residuo (1–118)", fontweight="bold")
    ax.set_ylabel("RMSF Cα (Å)", fontweight="bold")
    ax.set_ylim(0.0, ymax_rmsf)
    ax.set_title("B. Fluttuazione Atomica Locale RMSF(Cα) lungo la Sequenza (1–118)", fontweight="bold", fontsize=10.5)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper left", fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    f3_path = os.path.join(out_dir, f"{vname}_03_per_residue_qi_and_rmsf.pdf")
    plt.savefig(f3_path, dpi=300)
    plt.close(fig)

    # ==========================================================================
    # REPORT TESTUALE DELLA SINGOLA SIMULAZIONE: <VARIANT>_summary_report.txt
    # ==========================================================================
    rep_path = os.path.join(out_dir, f"{vname}_summary_report.txt")
    with open(rep_path, "w") as f:
        f.write("==========================================================================================\n")
        f.write(f" REPORT SINGOLA SIMULAZIONE: {title_label}\n")
        f.write("==========================================================================================\n")
        f.write(f"Numero totale di frame analizzati : {nf} (2.50 ns a 400 K)\n")
        f.write(f"Finestra equilibrata finale (30%) : Frame {f_start}..{nf-1} ({t_win_start:.2f}–2.50 ns)\n")
        f.write("------------------------------------------------------------------------------------------\n")
        f.write(" 1. METRICHE GLOBALI NELLA FINESTRA FINALE (30%):\n")
        f.write(f"    - Q: hydrophilic-all (Bekker r=0.84) : {np.mean(q_arr[f_start:,2]):.4f} ± {np.std(q_arr[f_start:,2]):.4f}\n")
        f.write(f"    - Q: all-all         (Bekker r=0.79) : {np.mean(q_arr[f_start:,1]):.4f} ± {np.std(q_arr[f_start:,1]):.4f}\n")
        f.write(f"    - Q: hydrophobic-small (r=0.76)      : {np.mean(q_arr[f_start:,3]):.4f} ± {np.std(q_arr[f_start:,3]):.4f}\n")
        f.write(f"    - Backbone RMSD Globale              : {np.mean(obs_arr[f_start:,1]):.2f} ± {np.std(obs_arr[f_start:,1]):.2f} Å\n")
        f.write(f"    - Functional Mask 3AAA RMSD          : {np.mean(obs_arr[f_start:,2]):.2f} ± {np.std(obs_arr[f_start:,2]):.2f} Å\n")
        f.write(f"    - Raggio di Girazione (Rgyr)         : {np.mean(obs_arr[f_start:,3]):.2f} ± {np.std(obs_arr[f_start:,3]):.2f} Å\n")
        f.write(f"    - SASA Totale                        : {np.mean(obs_arr[f_start:,4]):.1f} ± {np.std(obs_arr[f_start:,4]):.1f} Å²\n")
        f.write(f"    - % Alpha-Helix (STRIDE)             : {np.mean(obs_arr[f_start:,5]*100):.2f}%\n")
        if hb_arr is not None and len(hb_arr) == nf:
            f.write(f"    - H-Bonds Intra-Proteici Unici       : {np.mean(hb_arr[f_start:,1]):.1f} ± {np.std(hb_arr[f_start:,1]):.1f}\n")
        f.write("------------------------------------------------------------------------------------------\n")
        f.write(" 2. STABILITA' PER REGIONE STRUTTURALE (YANG ET AL. 1998) E INTERFACCE (FINESTRA FINALE 30%):\n")
        f.write(f"    - Q_Nterm (res 1–27)                 : {np.mean(q_arr[f_start:,4]):.4f}\n")
        f.write(f"    - Q_Rep1  (Ankyrin Rep 1, res 28–60) : {np.mean(q_arr[f_start:,5]):.4f}\n")
        f.write(f"    - Q_Rep2  (Ankyrin Rep 2, res 61–93) : {np.mean(q_arr[f_start:,6]):.4f}\n")
        f.write(f"    - Q_Cterm (res 94–118)               : {np.mean(q_arr[f_start:,7]):.4f}\n")
        f.write(f"    - Interfacce (N-R1 / R1-R2 / R2-C)   : {np.mean(q_arr[f_start:,8]):.3f} / {np.mean(q_arr[f_start:,9]):.3f} / {np.mean(q_arr[f_start:,10]):.3f}\n")
        f.write("------------------------------------------------------------------------------------------\n")
        f.write(f" 3. VALORI LOCALI SUI SITI CHIAVE (37, 60, 68, 102, 103):\n")
        for r_focus in [37, 60, 68, 102, 103]:
            idx = r_focus - 1
            f.write(f"    - Residuo {r_focus:3d} ({res_names[idx]:3s}): <Q_i>_400K = {qi_vals[idx]:.4f} | N_cont_0 = {n_cont[idx]:3d} | RMSF_CA = {rmsf_arr[idx,1]:.2f} Å\n")
        f.write("------------------------------------------------------------------------------------------\n")
        f.write(f" 4. RESIDUI DELLO SCAFFOLD LOCALMENTE INSTABILI/DEBOLI (N_0 >= 30):\n")
        f.write(f"    - Con <Q_i>_400K < 0.60 (Critici) : {n_sev_scaf} residui\n")
        f.write(f"    - Con <Q_i>_400K < 0.75 (Deboli)  : {n_unst_scaf} residui\n")
        for idx in np.where(unstable_scaf_mask)[0]:
            tag_s = "[CRITICO <0.60]" if qi_vals[idx] < 0.60 else "[DEBOLE  <0.75]"
            f.write(f"    * Res {res_ids[idx]:3d} ({res_names[idx]:3s}) -> <Q_i> = {qi_vals[idx]:.4f} {tag_s} (N_0 = {n_cont[idx]}, RMSF = {rmsf_arr[idx,1]:.2f} Å)\n")
        f.write("==========================================================================================\n")

    print(f"[OK] Plot e report individuali per {vname} salvati in: {out_dir}/")
    return True


if __name__ == "__main__":
    n_done = 0
    for spec in RUN_SPECS:
        if process_single_run(*spec):
            n_done += 1
    if n_done == 0:
        print("ATTENZIONE: Nessuna simulazione analizzata trovata in 05_analysis/. Esegui prima 09_analyze_matched_validation.tcl")
    else:
        print(f"\n>>> Generati con successo i pacchetti di plot individuali per {n_done} simulazioni in 05_analysis/single_runs/!")
