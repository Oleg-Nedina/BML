#!/usr/bin/env python3
"""
12_plot_hbonds_and_helices.py
Genera i grafici e il report numerico per le nuove metriche richieste dal Prof. Gautieri:
  1. Numero di Legami a Idrogeno (3.5 Å / 30°, totali proteina + backbone elicoidale)
     e occupazione (%) dei legami H chiave.
  2. Quantità di alfa-eliche e comportamento delle 7 singole eliche native (H1..H7),
     con focus speciale sulle ultime eliche (H6: 82-89, H7: 112-117 C-terminale)
     e sul modulo critico Repeat 1 (H3: 39-44, H4: 49-57).
  3. Supporta automaticamente sia n=1 (Seed 1) sia il triplicato (n=3 seeds -> μ ± σ!).

Salva i risultati in:
  - 06_presentation/presentazione_prof_WT_vs_Triplo/05_hbonds_comparison_WT_vs_Mutants.pdf (+ .png)
  - 06_presentation/presentazione_prof_WT_vs_Triplo/06_alpha_helices_timeline_and_H1_H7_behavior.pdf (+ .png)
  - 06_presentation/presentazione_prof_WT_vs_Triplo/05_hbonds_and_helices_summary.txt

Esecuzione dalla cartella BML/:
  python3 scripts/12_plot_hbonds_and_helices.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

IN_DIR = "05_analysis/helices_hbonds"
OUT_DIR = "06_presentation/presentazione_prof_WT_vs_Triplo"
os.makedirs(OUT_DIR, exist_ok=True)

VARIANT_ORDER = ["WT", "K37L", "T103L", "K37L_T103L", "K37L_H68G_T103L"]
VARIANT_LABELS = {
    "WT":              "WT (2MYO)",
    "K37L":            "K37L (Rep1)",
    "T103L":           "T103L (C-term)",
    "K37L_T103L":      "K37L+T103L (Doppio)",
    "K37L_H68G_T103L": "K37L+H68G+T103L (Triplo)",
}
VARIANT_COLORS = {
    "WT":              "#d62728",
    "K37L":            "#1f77b4",
    "T103L":           "#ff7f0e",
    "K37L_T103L":      "#2ca02c",
    "K37L_H68G_T103L": "#6a3d9a",
}

HELICES = [
    ("H1 (3–12)\nN-term",   0,  3, 12),
    ("H2 (15–23)\nN-term",  1, 15, 23),
    ("H3 (39–44)\nRep1",    2, 39, 44),
    ("H4 (49–57)\nRep1",    3, 49, 57),
    ("H5 (71–77)\nRep2",    4, 71, 77),
    ("H6 (82–89)\nRep2",    5, 82, 89),
    ("H7 (112–117)\nC-term",6,112,117),
]

data = {}
for vname in VARIANT_ORDER:
    gfiles = sorted(glob.glob(f"{IN_DIR}/{vname}_s*_helices_geom.dat"))
    if not gfiles:
        continue
    reps = []
    for gf in gfiles:
        pref = gf.replace("_helices_geom.dat", "")
        hb_p = np.loadtxt(f"{pref}_hbonds35_prot.dat")
        hb_b = np.loadtxt(f"{pref}_hbonds35_bb.dat")
        geom = np.loadtxt(gf)
        smat = np.loadtxt(f"{pref}_stride_matrix.dat")
        nf = len(geom)
        f_start = max(1, int(nf * 0.70))
        reps.append({
            "hb_p": hb_p,
            "hb_b": hb_b,
            "geom": geom,
            "smat": smat,
            "nf": nf,
            "f_start": f_start,
        })
    data[vname] = reps

if not data or "WT" not in data:
    print("ATTENZIONE: Nessun file trovato in 05_analysis/helices_hbonds/. Esegui prima: vmd -dispdev text -e scripts/12_analyze_hbonds_and_helices.tcl")
    raise SystemExit(0)

active_vars = [v for v in VARIANT_ORDER if v in data]


def save_both(fig, stem):
    pdf_p = os.path.join(OUT_DIR, f"{stem}.pdf")
    png_p = os.path.join(OUT_DIR, f"{stem}.png")
    fig.savefig(pdf_p, dpi=300, bbox_inches="tight")
    fig.savefig(png_p, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Salvato: {pdf_p} (+ .png)")


# ==============================================================================
# 1. REPORT NUMERICO: 05_hbonds_and_helices_summary.txt
# ==============================================================================
rep_path = os.path.join(OUT_DIR, "05_hbonds_and_helices_summary.txt")
with open(rep_path, "w") as f:
    f.write("========================================================================================================================\n")
    f.write(" ANALISI LEGAMI A IDROGENO (3.5 Å / 30°) E COMPORTAMENTO DELLE 7 ALFA-ELICHE NATIVE (H1–H7) A 400 K (FINESTRA FINALE 30%)\n")
    f.write(" Confini strutturali unificati: Yang et al. (1998) — N-term (1–27), Repeat 1 (28–60), Repeat 2 (61–93), C-term (94–118)\n")
    f.write("========================================================================================================================\n")
    f.write(f"{'Sistema':<18} {'Seeds':<6} {'H-Bonds Tot (3.5Å)':<20} {'H-Bonds BB (Eliche)':<20} {'H1%':<7} {'H2%':<7} {'H3%':<7} {'H4%':<7} {'H5%':<7} {'H6%':<7} {'H7%(C-cap)'}\n")
    f.write("-" * 120 + "\n")
    for vname in active_vars:
        reps = data[vname]
        hb_tot = [np.mean(r["hb_p"][r["f_start"]:, 1]) for r in reps]
        hb_bb  = [np.mean(r["hb_b"][r["f_start"]:, 1]) for r in reps]
        h_fracs = []
        for hi in range(7):
            col = 1 + hi * 3
            h_fracs.append(np.mean([np.mean(r["geom"][r["f_start"]:, col]) * 100.0 for r in reps]))
        f.write(
            f"{vname:<18s} {len(reps):<6d} "
            f"{np.mean(hb_tot):5.1f} ± {np.std(hb_tot) if len(reps)>1 else np.std(reps[0]['hb_p'][reps[0]['f_start']:,1]):4.1f}         "
            f"{np.mean(hb_bb):5.1f} ± {np.std(hb_bb) if len(reps)>1 else np.std(reps[0]['hb_b'][reps[0]['f_start']:,1]):4.1f}         "
            + " ".join(f"{hf:5.1f}%" for hf in h_fracs) + "\n"
        )

    f.write("\n------------------------------------------------------------------------------------------------------------------------\n")
    f.write(" DEFORMAZIONE (RMSD INTERNO CA, Å) E FRAYING (DISTANZA END-TO-END d_N-C, Å) DELLE ULTIME ELICHE (H6, H7) E DI REPEAT 1 (H3, H4):\n")
    f.write(f"{'Sistema':<18} {'H3 (39-44) RMSD / d_NC':<25} {'H4 (49-57) RMSD / d_NC':<25} {'H6 (82-89) RMSD / d_NC':<25} {'H7 (112-117) RMSD / d_NC'}\n")
    f.write("-" * 120 + "\n")
    for vname in active_vars:
        reps = data[vname]
        cols_s = []
        for hi in [2, 3, 5, 6]:  # H3, H4, H6, H7
            c_rmsd = 1 + hi * 3 + 1
            c_ee   = 1 + hi * 3 + 2
            rm = np.mean([np.mean(r["geom"][r["f_start"]:, c_rmsd]) for r in reps])
            ee = np.mean([np.mean(r["geom"][r["f_start"]:, c_ee]) for r in reps])
            cols_s.append(f"{rm:.2f} Å / {ee:.2f} Å")
        f.write(f"{vname:<18s} {cols_s[0]:<25s} {cols_s[1]:<25s} {cols_s[2]:<25s} {cols_s[3]}\n")
    f.write("========================================================================================================================\n")

print(f"[OK] Report H-Bonds & 7 Alfa-Eliche salvato in: {rep_path}")


# ==============================================================================
# 2. PLOT 1: 05_hbonds_comparison_WT_vs_Mutants
# ==============================================================================
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
fig.suptitle("Analisi della Rete di Legami a Idrogeno a 400 K (Cutoff Standard: 3.5 Å, Angolo 30°)",
             fontsize=13, fontweight="bold", y=0.99)

# A: H-Bonds Totali Intra-Proteici nel tempo
ax = axes[0]
ax.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85, label="Finestra finale 30%")
for vname in active_vars:
    reps = data[vname]
    t_ns = np.linspace(0.0, 2.50, reps[0]["nf"])
    y_m = np.mean([r["hb_p"][:, 1] for r in reps], axis=0)
    # Media mobile leggera (15 frame = 30 ps) per leggibilità
    k = np.ones(15) / 15.0
    y_sm = np.convolve(y_m, k, mode="same")
    lw = 2.2 if vname in ("WT", "K37L_H68G_T103L") else 1.4
    ax.plot(t_ns[8:-8], y_sm[8:-8], color=VARIANT_COLORS[vname], lw=lw, label=VARIANT_LABELS[vname])
ax.set_xlabel("Tempo a 400 K (ns)", fontweight="bold")
ax.set_ylabel("N° H-Bonds Totali (Proteina)", fontweight="bold")
ax.set_title("A. H-Bonds Totali nel Tempo (media mobile 30 ps)", fontweight="bold", fontsize=10)
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="lower left", fontsize=7.8)

# B: H-Bonds Backbone (Elicoidali) nel tempo
ax = axes[1]
ax.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85)
for vname in active_vars:
    reps = data[vname]
    t_ns = np.linspace(0.0, 2.50, reps[0]["nf"])
    y_m = np.mean([r["hb_b"][:, 1] for r in reps], axis=0)
    k = np.ones(15) / 15.0
    y_sm = np.convolve(y_m, k, mode="same")
    lw = 2.2 if vname in ("WT", "K37L_H68G_T103L") else 1.4
    ax.plot(t_ns[8:-8], y_sm[8:-8], color=VARIANT_COLORS[vname], lw=lw, label=VARIANT_LABELS[vname])
ax.set_xlabel("Tempo a 400 K (ns)", fontweight="bold")
ax.set_ylabel("N° H-Bonds Backbone (Eliche)", fontweight="bold")
ax.set_title("B. H-Bonds Backbone–Backbone nel Tempo", fontweight="bold", fontsize=10)
ax.grid(True, linestyle=":", alpha=0.5)

# C: Barplot comparativo (Finestra finale 30%)
ax = axes[2]
x = np.arange(len(active_vars))
w = 0.36
m_tot = [np.mean([np.mean(r["hb_p"][r["f_start"]:, 1]) for r in data[v]]) for v in active_vars]
m_bb  = [np.mean([np.mean(r["hb_b"][r["f_start"]:, 1]) for r in data[v]]) for v in active_vars]
b1 = ax.bar(x - w/2, m_tot, width=w, color="#4c72b0", edgecolor="black", label="H-Bonds Totali")
b2 = ax.bar(x + w/2, m_bb,  width=w, color="#55a868", edgecolor="black", label="H-Bonds Backbone")
for b, m in zip(b1, m_tot):
    ax.text(b.get_x()+b.get_width()/2, m+1.0, f"{m:.1f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold")
for b, m in zip(b2, m_bb):
    ax.text(b.get_x()+b.get_width()/2, m+1.0, f"{m:.1f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color="#1b611b")
ax.set_xticks(x)
ax.set_xticklabels([VARIANT_LABELS[v].replace(" (", "\n(") for v in active_vars], fontsize=7.8, fontweight="bold")
ax.set_ylabel("Media H-Bonds (Ultimi 30%)", fontweight="bold")
ax.set_ylim(0, max(m_tot) * 1.22)
ax.set_title("C. Media H-Bonds nella Finestra Finale 30%", fontweight="bold", fontsize=10)
ax.grid(True, axis="y", linestyle=":", alpha=0.5)
ax.legend(loc="upper right", fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.92])
save_both(fig, "05_hbonds_comparison_WT_vs_Mutants")


# ==============================================================================
# 3. PLOT 2: 06_alpha_helices_timeline_and_H1_H7_behavior
# ==============================================================================
fig = plt.figure(figsize=(15.0, 10.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.35, wspace=0.22)
fig.suptitle("Quantità di α-Eliche e Comportamento delle 7 Singole Eliche Native (H1–H7) a 400 K\n(Confini Strutturali Yang et al. 1998 — Focus su Repeat 1 [H3, H4] e Ultime Eliche [H6, H7])",
             fontsize=13, fontweight="bold", y=0.98)

cmap_ss = ListedColormap(["#f0f0f0", "#a6cee3", "#6a3d9a"])  # 0=Coil/Turn, 1=3_10, 2=Alpha-Helix

# A: Timeline STRIDE WT
ax1 = fig.add_subplot(gs[0, 0])
wt_mat = data["WT"][0]["smat"][:, 1:].T  # (118, n_t)
ax1.imshow(wt_mat, aspect="auto", origin="lower", cmap=cmap_ss, vmin=0, vmax=2, extent=[0, 2.5, 1, 118])
for b in (27.5, 60.5, 93.5):
    ax1.axhline(b, color="red", ls="--", lw=1.0)
for h_lbl, _, r1, r2 in HELICES:
    ax1.text(2.53, (r1+r2)/2.0, h_lbl.split("\n")[0], va="center", ha="left", fontsize=7.5, fontweight="bold")
ax1.set_title("A. Timeline STRIDE — Wild-Type (Viola = α-Elica, Grigio = Turn/Coil)", fontweight="bold", fontsize=10)
ax1.set_xlabel("Tempo a 400 K (ns)", fontweight="bold")
ax1.set_ylabel("Residuo (1–118)", fontweight="bold")

# B: Timeline STRIDE Triplo Mutante (se presente)
ax2 = fig.add_subplot(gs[0, 1])
tr_key = "K37L_H68G_T103L" if "K37L_H68G_T103L" in data else active_vars[-1]
tr_mat = data[tr_key][0]["smat"][:, 1:].T
ax2.imshow(tr_mat, aspect="auto", origin="lower", cmap=cmap_ss, vmin=0, vmax=2, extent=[0, 2.5, 1, 118])
for b in (27.5, 60.5, 93.5):
    ax2.axhline(b, color="red", ls="--", lw=1.0)
for h_lbl, _, r1, r2 in HELICES:
    ax2.text(2.53, (r1+r2)/2.0, h_lbl.split("\n")[0], va="center", ha="left", fontsize=7.5, fontweight="bold")
ax2.set_title(f"B. Timeline STRIDE — {VARIANT_LABELS[tr_key]} (Tenuta di H1–H7)", fontweight="bold", fontsize=10)
ax2.set_xlabel("Tempo a 400 K (ns)", fontweight="bold")
ax2.set_ylabel("Residuo (1–118)", fontweight="bold")

# C: Barplot Persistenza (%) delle 7 Eliche Native H1..H7 (Ultimi 30%)
ax3 = fig.add_subplot(gs[1, 0])
x = np.arange(7)
w_bar = 0.80 / len(active_vars)
for i, vname in enumerate(active_vars):
    reps = data[vname]
    vals = [np.mean([np.mean(r["geom"][r["f_start"]:, 1 + hi*3]) * 100.0 for r in reps]) for hi in range(7)]
    off = (i - (len(active_vars) - 1)/2.0) * w_bar
    ax3.bar(x + off, vals, width=w_bar*0.92, color=VARIANT_COLORS[vname], edgecolor="black", lw=0.4,
            label=VARIANT_LABELS[vname])
ax3.set_xticks(x)
ax3.set_xticklabels([h[0] for h in HELICES], fontsize=7.8, fontweight="bold")
ax3.set_ylabel("Persistenza α-Elica (%) — Ultimi 30%", fontweight="bold")
ax3.set_ylim(0, 112)
ax3.set_title("C. Conservazione di Ciascuna delle 7 α-Eliche Native (H1–H7)", fontweight="bold", fontsize=10)
ax3.grid(True, axis="y", linestyle=":", alpha=0.5)
ax3.legend(loc="lower left", fontsize=7.5, ncol=2)

# D: Comportamento delle Ultime Eliche (H6: 82-89 e H7: 112-117 C-cap) — Distanza End-to-End nel tempo
ax4 = fig.add_subplot(gs[1, 1])
ax4.axvspan(1.75, 2.50, color="#fff8dc", alpha=0.85, label="Finestra finale 30%")
for vname, ls_h6, ls_h7 in [("WT", "--", "-"), (tr_key, "--", "-")]:
    reps = data[vname]
    t_ns = np.linspace(0.0, 2.50, reps[0]["nf"])
    # H6 end-to-end (hi=5 -> col 1+5*3+2 = 18) & H7 end-to-end (hi=6 -> col 1+6*3+2 = 21)
    h6_ee = np.mean([r["geom"][:, 18] for r in reps], axis=0)
    h7_ee = np.mean([r["geom"][:, 21] for r in reps], axis=0)
    k = np.ones(15) / 15.0
    ax4.plot(t_ns[8:-8], np.convolve(h6_ee, k, mode="same")[8:-8], color=VARIANT_COLORS[vname], ls="--", lw=1.6,
             label=f"{vname} — H6 (82–89, media fin: {np.mean(h6_ee[reps[0]['f_start']:]):.1f} Å)")
    ax4.plot(t_ns[8:-8], np.convolve(h7_ee, k, mode="same")[8:-8], color=VARIANT_COLORS[vname], ls="-", lw=2.2,
             label=f"{vname} — H7 C-cap (112–117, fin: {np.mean(h7_ee[reps[0]['f_start']:]):.1f} Å)")
ax4.set_xlabel("Tempo a 400 K (ns)", fontweight="bold")
ax4.set_ylabel("Distanza End-to-End Elica d(Cα_N–Cα_C) (Å)", fontweight="bold")
ax4.set_title("D. Comportamento (Fraying / Stretching) delle Ultime Eliche H6 e H7", fontweight="bold", fontsize=10)
ax4.grid(True, linestyle=":", alpha=0.5)
ax4.legend(loc="best", fontsize=7.8)

save_both(fig, "06_alpha_helices_timeline_and_H1_H7_behavior")
print(f"\n>>> ANALISI H-BONDS E 7 ALFA-ELICHE COMPLETATA! Controlla: {OUT_DIR}/")
