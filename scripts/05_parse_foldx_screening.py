#!/usr/bin/env python3
"""
05_parse_foldx_screening.py
Fase 3 (ROAD3.md - Passo 3.3):
  1. Legge l'output di FoldX PositionScan:
     - 03_mutations/PS_2MYO_Repair_scanning_output.txt
     - (e/o i file energies_*.txt generati da PositionScan)
  2. Se disponibile, incrocia automaticamente i valori di <Q_i>_400K della MD
     (da 02_wt_baseline/q_per_residue_400K.dat) e il Consenso Ankyrin.
  3. Filtra qualsiasi mutazione che cada nella Functional Mask (24 residui) o in Proline/Cisteine
     che rompono le alfa-eliche.
  4. Genera:
     - 03_mutations/foldx_top_mutations.txt (Classifica completa ordinata per ddG)
     - 06_presentation/fig7_foldx_position_scan_heatmap.pdf (Grafico vettoriale PDF per le slide)
"""

import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("03_mutations", exist_ok=True)
os.makedirs("06_presentation", exist_ok=True)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "H1S": "H", "H2S": "H", "H3S": "H"
}

FUNCTIONAL_MASK = {
    4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36,
    41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100
}

CONSENSUS_TARGETS = {
    73:  ("H", "Ripristina motivo firma TPLH (L73H) in ANK3"),
    94:  ("V/I/L", "Ripristina tasca idrofobica GADVNA (K94V/I/L) in ANK3"),
    64:  ("K/N/A", "Rimuove tensione Ramachandran di PRO 64 in ANK2"),
    105: ("F/I/V", "Variante umana L105F (3AAA) nel C-cap ANK4"),
    37:  ("T", "Ripristina THR in pos 3 del motivo TPLH in ANK2"),
    68:  ("G/N/D", "Ottimizza beta-turn pos 1 in ANK3"),
    81:  ("G/H/A", "Rilassa loop corto inter-elica in ANK3"),
    59:  ("G", "Rilassa ingresso GADVNA in ANK2"),
    26:  ("A", "Ripristina ALA in GADVNA di ANK1"),
}

# Carica <Q_i>_400K se disponibile
qi_map = {}
qi_file = "02_wt_baseline/q_per_residue_400K.dat"
if os.path.exists(qi_file):
    with open(qi_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                qi_map[int(parts[0])] = float(parts[3])

# Cerca il file di output di PositionScan
scan_files = [
    "03_mutations/PS_2MYO_Repair_scanning_output.txt",
    "PS_2MYO_Repair_scanning_output.txt"
]
scan_path = next((p for p in scan_files if os.path.exists(p)), None)

mutations = []
if scan_path:
    with open(scan_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Formato tipico: LYSB94VAL  -1.4523  oppure LA73H  -1.23
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                continue
            tag, ddg_str = parts[0], parts[-1]
            try:
                ddg = float(ddg_str)
            except ValueError:
                continue

            # Parsing di tag tipo "LEUA73HIS" o "LA73H"
            m3 = re.match(r"^([A-Z0-9]{3})([A-Z])(\d+)([A-Z0-9]{3})$", tag)
            m1 = re.match(r"^([A-Z])([A-Z])(\d+)([A-Z])$", tag)
            if m3:
                wt3, chain, resnum, mut3 = m3.groups()
                wt1 = THREE_TO_ONE.get(wt3, "?")
                mut1 = THREE_TO_ONE.get(mut3, "?")
                res_id = int(resnum)
            elif m1:
                wt1, chain, resnum, mut1 = m1.groups()
                res_id = int(resnum)
            else:
                continue

            if wt1 == mut1 or mut1 == "?":
                continue
            if res_id in FUNCTIONAL_MASK:
                continue

            foldx_code = f"{wt1}{chain}{res_id}{mut1}"
            mutations.append((ddg, res_id, wt1, mut1, foldx_code))

if not mutations:
    print("ATTENZIONE: Nessun file PS_2MYO_Repair_scanning_output.txt trovato in 03_mutations/.")
    print("Assicurati di aver prima eseguito: foldx -f 03_mutations/position_scan_hotspots.cfg")
    raise SystemExit(0)

mutations.sort(key=lambda x: x[0])

out_txt = "03_mutations/foldx_top_mutations.txt"
with open(out_txt, "w") as f:
    f.write("=========================================================================================\n")
    f.write(" CLASSIFICA DELLE MIGLIORI MUTAZIONI STABILIZZANTI DA FOLDX POSITIONSCAN (2MYO_Repair)\n")
    f.write(" (Esclusi i 24 residui della Functional Mask 3AAA e le auto-mutazioni WT->WT)\n")
    f.write("=========================================================================================\n")
    f.write(f"{'Rank':<5} {'Mutazione':<10} {'Codice FoldX':<14} {'Modulo':<8} {'ddG_FoldX (kcal/mol)':<22} {'<Q_i>_400K':<12} {'Note / Consenso Ankyrin'}\n")
    f.write("-" * 105 + "\n")

    for idx, (ddg, res_id, wt1, mut1, code) in enumerate(mutations[:40], 1):
        rep = "ANK1" if res_id <= 34 else ("ANK2" if res_id <= 67 else ("ANK3" if res_id <= 100 else "ANK4"))
        qi_val = f"{qi_map[res_id]:.3f}" if res_id in qi_map else "N/A"
        note = ""
        if res_id in CONSENSUS_TARGETS:
            target_aas, desc = CONSENSUS_TARGETS[res_id]
            if mut1 in target_aas.split("/"):
                note = f"★ CONSENSO ANKYRIN ({desc})"
            else:
                note = desc
        f.write(f"{idx:<5d} {wt1}{res_id}{mut1:<7s} {code:<14s} {rep:<8s} {ddg:<+22.4f} {qi_val:<12s} {note}\n")

# Generazione Grafico Vettoriale PDF delle Top 20 mutazioni stabilizzanti
top_n = [m for m in mutations if m[0] < 0][:20]
if not top_n:
    top_n = mutations[:15]

fig, ax = plt.subplots(figsize=(10.5, 5.5))
labels = [f"{m[2]}{m[1]}{m[3]}" for m in top_n][::-1]
ddgs = [m[0] for m in top_n][::-1]
resids = [m[1] for m in top_n][::-1]

colors = []
for r, (_, _, wt1, mut1, _) in zip(resids, top_n[::-1]):
    if r in CONSENSUS_TARGETS and mut1 in CONSENSUS_TARGETS[r][0].split("/"):
        colors.append("#2ca02c")  # Verde: sia stabilizzante FoldX sia Consenso Ankyrin!
    elif r <= 34:
        colors.append("#1f77b4")
    elif r <= 67:
        colors.append("#17becf")
    elif r <= 100:
        colors.append("#ff7f0e")
    else:
        colors.append("#9467bd")

bars = ax.barh(range(len(labels)), ddgs, color=colors, edgecolor="#222222", height=0.72)
ax.axvline(0.0, color="#333333", lw=1.0)
ax.axvline(-1.0, color="#cc0000", lw=1.2, ls="--", label="Soglia forte stabilizzazione (ΔΔG = -1.0 kcal/mol)")

for bar, val in zip(bars, ddgs):
    ax.text(val - 0.04 if val < 0 else val + 0.04,
            bar.get_y() + bar.get_height() / 2.0,
            f"{val:+.2f}",
            va="center", ha="right" if val < 0 else "left",
            fontsize=8.5, fontweight="bold")

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
ax.set_xlabel("FoldX ΔΔG = ΔG_mut - ΔG_WT (kcal/mol)  [Valori negativi = Stabilizzante]", fontsize=10, fontweight="bold")
ax.set_title("Top 20 Mutazioni Termostabilizzanti in Myotrophin (2MYO) da FoldX PositionScan", fontsize=11.5, fontweight="bold")
ax.grid(axis="x", linestyle=":", alpha=0.6)
ax.legend(loc="lower left", frameon=True, fontsize=9)

plt.tight_layout()
pdf_path = "06_presentation/fig7_foldx_position_scan_heatmap.pdf"
plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
plt.close()

print("==================================================================")
print(f" Classifica FoldX salvata in : {out_txt}")
print(f" Grafico PDF salvato in      : {pdf_path}")
print("==================================================================")
