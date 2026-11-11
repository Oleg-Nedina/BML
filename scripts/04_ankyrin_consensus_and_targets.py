#!/usr/bin/env python3
"""
04_ankyrin_consensus_and_targets.py
Fase 3 (ROAD3.md - Passo 3.1):
  1. Allinea i 4 Ankyrin Repeats di Myotrophin (2MYO, 118 aa) contro il
     Consenso Universale Termostabile degli Ankyrin Repeats (Mosavi et al., PNAS 2002).
  2. Estrae i residui con angoli di Ramachandran anomali dal REMARK 500 di 2MYO.pdb.
  3. Applica la Functional Mask oggettiva di 24 residui (contatto < 4.5 A in 3AAA.pdb).
  4. Genera:
     - 03_mutations/ankyrin_consensus_report.txt (Report testuale completo)
     - 03_mutations/position_scan_hotspots.cfg   (Configurazione FoldX PositionScan)
     - 06_presentation/fig6_ankyrin_consensus_alignment.pdf (Grafico vettoriale PDF per le slide)
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

os.makedirs("03_mutations", exist_ok=True)
os.makedirs("06_presentation", exist_ok=True)

# 1. Sequenza esatta di 2MYO (118 aa, verificata da SEQRES e getFasta in Fase 1)
SEQ_2MYO = (
    "MCDKEFMWALKNGDLDEVKDYVAKGEDVNRTLE"  # 1..33  (ANK1)
    "GGRKPLHYAADCGQLEILEFLLLKGADINAPDK"  # 34..66 (ANK2)
    "HHITPLLSAVYEGHVSCVKLLLSKGADKTVKGP"  # 67..99 (ANK3)
    "DGLTAFEATDNQAIKALLQ"                # 100..118 (ANK4 - C-cap, nota: 105 è L in 2MYO)
)
# Correzione indice 104 (residuo 105: LEU in 2MYO ratto, PHE in 3AAA umano)
SEQ_2MYO = (
    "MCDKEFMWALKNGDLDEVKDYVAKGEDVNRTLE"  # 1..33
    "GGRKPLHYAADCGQLEILEFLLLKGADINAPDK"  # 34..66
    "HHITPLLSAVYEGHVSCVKLLLSKGADKTVKGP"  # 67..99
    "DGLTALEATDNQAIKALLQ"                # 100..118 (res 105 = L)
)
assert len(SEQ_2MYO) == 118, f"Errore lunghezza sequenza: {len(SEQ_2MYO)}"

# 2. I 24 residui della Functional Mask (< 4.5 A da Actin-Capping Protein in 3AAA)
FUNCTIONAL_MASK = {
    4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36,
    41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100
}

# 3. Residui con angoli di Ramachandran fuori regione in 2MYO.pdb (REMARK 500, righe 126-160)
RAMA_OUTLIERS = {
    11, 24, 26, 27, 28, 29, 31, 32, 33, 38, 45, 47, 48, 53, 54, 57,
    59, 60, 62, 64, 68, 77, 81, 84, 85, 87, 88, 92, 93, 94, 97, 100,
    103, 105, 108
}

# 4. Consenso Universale Ankyrin Repeat di 33 aa (Mosavi et al., PNAS 2002; Binz et al. DARPins)
# Posizioni 1..33 del mattoncino canonico:
# Posizione 1 nel repeat corrisponde a:
#   - ANK1: resid 2..34  (con gap iniziale sul loop N-cap)
#   - ANK2: resid 35..67 (G35=pos1, R36=pos2, K37=pos3, P38=pos4, L39=pos5, H40=pos6...)
#   - ANK3: resid 68..100 (H68=pos1, I69=pos2, T70=pos3, P71=pos4, L72=pos5, L73=pos6...)
#   - ANK4: resid 101..118 (G101=pos1, L102=pos2, T103=pos3, A104=pos4, L105=pos5...)
ANK_CONSENSUS = "GXTPLHLAAXXGHLEIVEVLLKXGADVNAXDXX"

# Mappatura dei 4 repeat di 2MYO sulle 33 colonne del consenso canonico (resid 1-based, o None se gap)
REPEAT_MAP = {
    "ANK1 (N-cap, 2-34)":  [None, None, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34][:33],
    "ANK2 (Core, 35-67)":  list(range(35, 68)),
    "ANK3 (Core, 68-100)": list(range(68, 101)),
    "ANK4 (C-cap, 101-118)": list(range(101, 119)) + [None] * 15,
}
# Allineamento esatto strutturale per ANK1 (Elica 1A: 5-11, Loop: 12-14, Elica 1B: 15-22, Loop GADVNA: 23-34)
REPEAT_MAP["ANK1 (N-cap, 2-34)"] = [
    None, None, None, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, None, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32
]

# Raccogliamo le deviazioni chiave dal consenso (fuori dalla Functional Mask)
KEY_HOTSPOTS = [
    (73,  "L", "H", "ANK3 (pos 6)",  "Motivo firma TPLH interrotto (T70-P71-L72-L73): L73H ripristina l'Istidina canonica di consenso!"),
    (94,  "K", "V/I", "ANK3 (pos 27)", "Motivo GADVNA interrotto (G91-A92-D93-K94): Lisina carica in tasca idrofobica V/I + Ramachandran outlier!"),
    (64,  "P", "K/N/A", "ANK2 (pos 30)", "Loop GADINAPD: PRO 64 ha il peggior Ramachandran outlier di 2MYO (psi=+121.5, phi=-27.7)!"),
    (105, "L", "F", "ANK4 (pos 5)",  "Unica differenza evolutiva Ratto (LEU 105) vs Uomo (PHE 105 in 3AAA) + Ramachandran outlier!"),
    (37,  "K", "T", "ANK2 (pos 3)",  "Posizione 3 del motivo TPLH occupata da LYS 37 (invece di THR consenso; G35-R36-K37-P38-L39-H40)"),
    (68,  "H", "G/N", "ANK3 (pos 1)", "Posizione 1 beta-turn (di solito GLY) occupata da HIS 68 con Ramachandran positivo anomalo"),
    (81,  "V", "G/H", "ANK3 (pos 14)", "Loop corto inter-elica in ANK3: VAL 81 ha angoli di Ramachandran anomali (psi=+4.0, phi=+38.4)"),
    (59,  "A", "G", "ANK2 (pos 25)", "Ingresso del loop GADVNA: ALA 59 e' in tensione sterica rispetto a GLY/ALA"),
    (26,  "E", "A", "ANK1 (pos 25)", "Motivo GADVNA in ANK1 (G23-K24-G25-E26-D27-V28-N29): GLU 26 ha Ramachandran outlier"),
    (103, "T", "T/P", "ANK4 (pos 3)", "Ingresso elica C-cap ANK4 (G101-L102-T103-A104-L105): THR 103 Ramachandran outlier"),
]

# Scrittura del report dettagliato
report_path = "03_mutations/ankyrin_consensus_report.txt"
with open(report_path, "w") as f:
    f.write("====================================================================================\n")
    f.write(" REPORT FASE 3.1: ALLINEAMENTO CONSENSO ANKYRIN REPEATS & TARGET PER FOLDX (2MYO)\n")
    f.write(" Riferimento: Mosavi, Minor & Peng (PNAS 2002, 99:16029-16034)\n")
    f.write("====================================================================================\n\n")

    f.write("1. ALLINEAMENTO DEI 4 MODULI ANKYRIN DI 2MYO CONTRO IL CONSENSO (33 POSIZIONI)\n")
    f.write("   Legenda: [x] = Residuo nella Functional Mask (3AAA < 4.5 A, NON mutabile)\n")
    f.write("            *   = Angolo di Ramachandran fuori regione in 2MYO.pdb (REMARK 500)\n\n")

    pos_header = "Pos Repeat : " + " ".join(f"{i:3d}" for i in range(1, 34))
    con_header = "CONSENSO   : " + " ".join(f"  {c}" for c in ANK_CONSENSUS)
    f.write(pos_header + "\n")
    f.write(con_header + "\n")
    f.write("-" * len(pos_header) + "\n")

    for rname, rlist in REPEAT_MAP.items():
        row_aa = []
        for res_id in rlist:
            if res_id is None:
                row_aa.append("  -")
            else:
                aa = SEQ_2MYO[res_id - 1]
                if res_id in FUNCTIONAL_MASK:
                    row_aa.append(f"[{aa}]")
                elif res_id in RAMA_OUTLIERS:
                    row_aa.append(f" *{aa}")
                else:
                    row_aa.append(f"  {aa}")
        f.write(f"{rname:22s}: " + " ".join(row_aa) + "\n")

    f.write("\n====================================================================================\n")
    f.write("2. TOP 10 HOTSPOTS STRUTTURALI & EVOLUTIVI (FUORI DALLA FUNCTIONAL MASK)\n")
    f.write("====================================================================================\n")
    f.write(f"{'Res':<6} {'WT':<4} {'Consenso':<10} {'Modulo':<16} {'Razionale Biofisico / Strutturale'}\n")
    f.write("-" * 100 + "\n")
    for res_id, wt_aa, target_aa, mod, reason in KEY_HOTSPOTS:
        assert res_id not in FUNCTIONAL_MASK, f"ERRORE: {res_id} e' nella Functional Mask!"
        assert SEQ_2MYO[res_id - 1] == wt_aa, f"Discrepanza residuo {res_id}: {SEQ_2MYO[res_id-1]} vs {wt_aa}"
        f.write(f"{res_id:<6d} {wt_aa:<4s} {target_aa:<10s} {mod:<16s} {reason}\n")

    # Costruiamo la lista di posizioni per FoldX PositionScan:
    # Uniamo i 10 Hotspot di consenso + tutti i Ramachandran outliers non-mascherati e non-terminali (resid 4..115)
    scan_resids = sorted(list(
        ({h[0] for h in KEY_HOTSPOTS} | RAMA_OUTLIERS | {35, 40, 70, 71, 72, 91, 95, 96, 101, 104, 110})
        - FUNCTIONAL_MASK
        - {1, 2, 3, 116, 117, 118}
    ))

    foldx_codes = [f"{SEQ_2MYO[r-1]}A{r}d" for r in scan_resids]
    f.write("\n====================================================================================\n")
    f.write(f"3. LISTA DI {len(scan_resids)} POSIZIONI CANDIDATE PER FOLDX POSITIONSCAN (20 AA ciascuna)\n")
    f.write("====================================================================================\n")
    f.write("Stringa FoldX positions:\n")
    f.write(",".join(foldx_codes) + "\n")

# Scrittura del file di configurazione FoldX PositionScan
cfg_path = "03_mutations/position_scan_hotspots.cfg"
with open(cfg_path, "w") as f:
    f.write("command=PositionScan\n")
    f.write("pdb-dir=03_mutations\n")
    f.write("output-dir=03_mutations\n")
    f.write("pdb=2MYO_Repair.pdb\n")
    f.write("positions=" + ",".join(foldx_codes) + "\n")

# Generazione Grafico Vettoriale PDF (06_presentation/fig6_ankyrin_consensus_alignment.pdf)
fig, ax = plt.subplots(figsize=(15, 4.8))

rows_labels = [
    "CONSENSO ANK\n(Mosavi et al. 2002)",
    "ANK1 (N-cap, 2-34)",
    "ANK2 (Core, 35-67)",
    "ANK3 (Core, 68-100)",
    "ANK4 (C-cap, 101-118)"
]

n_cols = 33
n_rows = 5

for r_idx in range(n_rows):
    y = n_rows - 1 - r_idx
    if r_idx == 0:
        for c_idx, c_aa in enumerate(ANK_CONSENSUS):
            color = "#2b5c8f" if c_aa != "X" else "#d9e1e8"
            tcolor = "white" if c_aa != "X" else "#555555"
            rect = mpatches.FancyBboxPatch((c_idx + 0.08, y + 0.12), 0.84, 0.76,
                                           boxstyle="round,pad=0.05", facecolor=color, edgecolor="none")
            ax.add_patch(rect)
            ax.text(c_idx + 0.5, y + 0.5, c_aa, ha="center", va="center",
                    color=tcolor, fontweight="bold", fontsize=9.5)
    else:
        rkey = list(REPEAT_MAP.keys())[r_idx - 1]
        rlist = REPEAT_MAP[rkey]
        for c_idx, res_id in enumerate(rlist):
            if res_id is None:
                rect = mpatches.FancyBboxPatch((c_idx + 0.08, y + 0.12), 0.84, 0.76,
                                               boxstyle="round,pad=0.05", facecolor="#f2f2f2", edgecolor="#dddddd")
                ax.add_patch(rect)
                ax.text(c_idx + 0.5, y + 0.5, "-", ha="center", va="center", color="#999999", fontsize=9)
                continue

            aa = SEQ_2MYO[res_id - 1]
            cons_aa = ANK_CONSENSUS[c_idx]
            if res_id in FUNCTIONAL_MASK:
                fcolor = "#ffb3b3"   # Rosso chiaro: Functional Mask (intoccabile)
                ecolor = "#cc0000"
            elif res_id in {h[0] for h in KEY_HOTSPOTS}:
                fcolor = "#ffd966"   # Giallo/Oro: Hotspot primario (deviante dal consenso / Ramachandran)
                ecolor = "#d68a00"
            elif cons_aa != "X" and aa == cons_aa:
                fcolor = "#b6d7a8"   # Verde: Identico al consenso Ankyrin
                ecolor = "#38761d"
            else:
                fcolor = "#ffffff"
                ecolor = "#bbbbbb"

            rect = mpatches.FancyBboxPatch((c_idx + 0.08, y + 0.12), 0.84, 0.76,
                                           boxstyle="round,pad=0.05", facecolor=fcolor, edgecolor=ecolor, linewidth=1.1)
            ax.add_patch(rect)
            ax.text(c_idx + 0.5, y + 0.60, aa, ha="center", va="center",
                    color="#111111", fontweight="bold", fontsize=9)
            ax.text(c_idx + 0.5, y + 0.27, f"{res_id}", ha="center", va="center",
                    color="#444444", fontsize=6.5)

# Annotazioni delle regioni strutturali canoniche sopra il grafico
regions = [
    (0, 4, "Beta-Turn", "#e6f2ff"),
    (4, 12, "Helix A (Inner, TPLH)", "#d0e8ff"),
    (12, 14, "Loop", "#f0f0f0"),
    (14, 23, "Helix B (Outer)", "#d0e8ff"),
    (23, 33, "Inter-Repeat Loop (GADVNA)", "#fff2cc")
]
for c_start, c_end, label, rcol in regions:
    ax.annotate(label, xy=((c_start + c_end) / 2.0, n_rows + 0.15),
                ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#222222")
    ax.plot([c_start + 0.1, c_end - 0.1], [n_rows + 0.08, n_rows + 0.08], color="#444444", lw=1.5)

ax.set_xlim(0, n_cols)
ax.set_ylim(-0.1, n_rows + 0.55)
ax.set_xticks([i + 0.5 for i in range(n_cols)])
ax.set_xticklabels([str(i + 1) for i in range(n_cols)], fontsize=8)
ax.set_yticks([n_rows - 1 - i + 0.5 for i in range(n_rows)])
ax.set_yticklabels(rows_labels, fontsize=9, fontweight="bold")
ax.set_xlabel("Posizione Canonica all'interno del Modulo Ankyrin Repeat (1..33)", fontsize=10, fontweight="bold")
ax.set_title("Allineamento Strutturale dei 4 Ankyrin Repeats di Myotrophin (2MYO) vs Consenso Universale",
             fontsize=11.5, fontweight="bold", pad=18)

for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(axis="both", length=0)

legend_patches = [
    mpatches.Patch(facecolor="#b6d7a8", edgecolor="#38761d", label="Già identico al Consenso Ankyrin"),
    mpatches.Patch(facecolor="#ffd966", edgecolor="#d68a00", label="Hotspot Candidato (Deviazione dal Consenso / Tensione Ramachandran)"),
    mpatches.Patch(facecolor="#ffb3b3", edgecolor="#cc0000", label="Functional Mask (< 4.5 Å da Capping Protein in 3AAA — NON mutabile)"),
    mpatches.Patch(facecolor="#ffffff", edgecolor="#bbbbbb", label="Posizione variabile (X) o neutra dello scaffold"),
]
ax.legend(handles=legend_patches, loc="upper center", bbox_to_anchor=(0.5, -0.15),
          ncol=2, frameon=True, fontsize=8.5)

plt.tight_layout()
pdf_out = "06_presentation/fig6_ankyrin_consensus_alignment.pdf"
plt.savefig(pdf_out, format="pdf", bbox_inches="tight")
plt.close()

print("==================================================================")
print(" Analisi Consenso Ankyrin & Target FoldX completata con successo!")
print(f"  -> Report testuale : {report_path}")
print(f"  -> Config FoldX    : {cfg_path} ({len(scan_resids)} posizioni x 20 aa = {len(scan_resids)*20} mutazioni)")
print(f"  -> Grafico PDF     : {pdf_out}")
print("==================================================================")
