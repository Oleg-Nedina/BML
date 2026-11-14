#!/usr/bin/env python3
"""
06_build_consensus_matrix.py
Fase 3 (ROAD3.md - Passo 3.4 & 3.5):
Costruisce la Matrice di Consenso Multi-Tool ("Orthogonal Consensus Funnel")
che fonde i 5 Pilastri di un lavoro di ricerca:
  1. Dinamica Molecolare WT a 400 K (<Q_i>_400K e RMSF)
  2. FoldX 5 PositionScan + Scomposizione Energetica (energies_*_2MYO_Repair.txt:
     HBond_Backbone, HBond_Sidechain, VdW, Electrostatics, Solv_Hydrophobic, Clash, Torsion)
  3. Consenso Evolutivo degli Ankyrin Repeats (Mosavi et al., PNAS 2002) & Ratto/Uomo (3AAA)
  4. ThermoMPNN (se presente 03_mutations/thermompnn_2MYO.csv)
  5. FireProt 2.0 / Rosetta & DynaMut2 (se presenti 03_mutations/fireprot_results.csv / dynamut2_results.csv)

Genera:
  - 03_mutations/consensus_master_table.txt (Tabella Master completa + Top 3 Candidati suggeriti)
  - 06_presentation/fig8_consensus_energy_breakdown.pdf (Grafico vettoriale PDF della scomposizione energetica)
"""

import os
import re
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("03_mutations", exist_ok=True)
os.makedirs("06_presentation", exist_ok=True)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

FUNCTIONAL_MASK = {
    4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36,
    41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100
}

# Residui in alfa-elica interna/esterna (HELIX record di 2MYO.pdb) dove PRO/GLY rompono l'elica
HELIX_CORE_RESIDS = (
    set(range(4, 11)) | set(range(16, 22)) |
    set(range(40, 46)) | set(range(50, 56)) |
    set(range(72, 78)) | set(range(83, 88))
)

CONSENSUS_INFO = {
    (59, "G"):  ("★ Consenso ANK (GADVNA)", +2.0),
    (68, "G"):  ("★ Consenso ANK (Beta-Turn pos 1)", +2.0),
    (68, "N"):  ("★ Consenso ANK (Beta-Turn pos 1)", +1.8),
    (81, "L"):  ("★ Packing Loop Inter-Elica ANK3", +1.5),
    (103, "F"): ("★ Core Hydrophobic C-cap ANK4", +1.8),
    (103, "L"): ("★ Consenso ANK (TPLH pos 5)", +1.8),
    (108, "L"): ("★ Packing Elica C-cap ANK4", +1.5),
    (105, "F"): ("★ Variante Umana Naturale (3AAA)", +1.8),
    (37, "L"):  ("★ Stabilizzazione Motivo TPLH ANK2", +1.2),
    (37, "I"):  ("★ Stabilizzazione Motivo TPLH ANK2", +1.2),
    (60, "L"):  ("★ Rilascio Tensione Loop GADVNA ANK2", +1.2),
    (64, "V"):  ("★ Rimuove Outlier Ramachandran PRO 64", +1.5),
    (93, "G"):  ("★ Rilassa Loop GADKTV in ANK3", +1.2),
    (27, "F"):  ("★ Packing Loop GADVNA in ANK1", +1.0),
}

# 1. Carica <Q_i>_400K e RMSF_400K dalla MD del Wild-Type
qi_map = {}
if os.path.exists("02_wt_baseline/q_per_residue_400K.dat"):
    with open("02_wt_baseline/q_per_residue_400K.dat") as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                p = line.split()
                qi_map[int(p[0])] = float(p[3])

rmsf_map = {}
if os.path.exists("02_wt_baseline/rmsf_400K.dat"):
    with open("02_wt_baseline/rmsf_400K.dat") as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                p = line.split()
                rmsf_map[int(p[0])] = float(p[1])

# 2. Carica ThermoMPNN / FireProt / DynaMut2 se l'utente ha scaricato i CSV in 03_mutations/
mpnn_map = {}
if os.path.exists("03_mutations/thermompnn_2MYO.csv"):
    with open("03_mutations/thermompnn_2MYO.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Supporta colonne standard ThermoMPNN (Mutation / mut / ddG_pred / ddG)
            mut_str = (row.get("Mutation") or row.get("mutation") or row.get("mut") or "").strip()
            ddg_val = row.get("ddG_pred") or row.get("ddG") or row.get("ddg")
            m = re.search(r"([A-Z])(\d+)([A-Z])$", mut_str)
            if m and ddg_val is not None:
                try:
                    mpnn_map[(int(m.group(2)), m.group(3))] = float(ddg_val)
                except ValueError:
                    pass

fireprot_map = {}
if os.path.exists("03_mutations/fireprot_results.csv"):
    with open("03_mutations/fireprot_results.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Supporta formato tabellare FireProt (colonna mutation o chain+position+mutant, ddG Rosetta/FoldX)
            mut_str = (row.get("Mutation") or row.get("mutation") or "").strip()
            pos_str = (row.get("Position") or row.get("position") or row.get("res_num") or "").strip()
            sub_str = (row.get("Substitution") or row.get("mutant") or row.get("mut_aa") or "").strip()
            ddg_val = (
                row.get("Rosetta") or row.get("rosetta_ddg") or row.get("ddG_Rosetta")
                or row.get("FoldX") or row.get("foldx_ddg") or row.get("ddG")
            )
            res_num, mut_aa = None, None
            m = re.search(r"([A-Z])(\d+)([A-Z])", mut_str)
            if m:
                res_num, mut_aa = int(m.group(2)), m.group(3)
            elif pos_str.isdigit() and len(sub_str) == 1:
                res_num, mut_aa = int(pos_str), sub_str.upper()
            if res_num is not None and mut_aa is not None:
                try:
                    fireprot_map[(res_num, mut_aa)] = float(ddg_val) if ddg_val not in (None, "", "-") else -1.0
                except ValueError:
                    fireprot_map[(res_num, mut_aa)] = -1.0

dynamut_map = {}
if os.path.exists("03_mutations/dynamut2_results.csv"):
    with open("03_mutations/dynamut2_results.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mut_str = (row.get("mutation") or row.get("Mutation") or "").strip()
            ddg_val = row.get("prediction") or row.get("ddG") or row.get("ddg")
            m = re.search(r"([A-Z])(\d+)([A-Z])", mut_str)
            if m and ddg_val is not None:
                try:
                    # In DynaMut2 valori POSITIVI sono stabilizzanti -> invertiamo in convenzione ΔΔG < 0 = stabilizzante
                    dynamut_map[(int(m.group(2)), m.group(3))] = -float(ddg_val)
                except ValueError:
                    pass

# 3. Estrai scomposizione energetica dettagliata da 03_mutations/energies_<res>_2MYO_Repair.txt
# Colonne FoldX energies_*.txt:
# 0:Label, 1:Total, 2:BackHbond, 3:SideHbond, 4:VanDerWaals, 5:Electrostatics,
# 6:SolvationPolar, 7:SolvationHydrophobic, 8:VdWClash, 9:EntropySidec, 10:EntropyMainc,
# 14:TorsionClash, 15:BackboneClash
energy_terms = {}
for efile in glob.glob("03_mutations/energies_*_2MYO_Repair.txt"):
    m_res = re.search(r"energies_(\d+)_2MYO_Repair\.txt$", efile)
    if not m_res:
        continue
    res_id = int(m_res.group(1))
    with open(efile) as f:
        lines = [l.strip().split("\t") for l in f if l.strip()]
    if not lines or not lines[0][0].startswith("WTref"):
        continue
    wt_vals = [float(x) for x in lines[0][1:]]
    for row in lines[1:]:
        aa3 = row[0].split("_")[0]
        if aa3 not in THREE_TO_ONE:
            continue
        aa1 = THREE_TO_ONE[aa3]
        vals = [float(x) for x in row[1:]]
        diff = [vals[i] - wt_vals[i] for i in range(len(wt_vals))]
        energy_terms[(res_id, aa1)] = {
            "dTotal": diff[0],
            "dHBond": diff[1] + diff[2],
            "dVdW": diff[3],
            "dElec": diff[4],
            "dSolvP": diff[5],
            "dSolvH": diff[6],
            "dClash": diff[7] + diff[13] + diff[14],
            "dEntropy": diff[8] + diff[9],
        }

# 4. Leggi le mutazioni da PS_2MYO_Repair_scanning_output.txt
records = []
scan_path = "03_mutations/PS_2MYO_Repair_scanning_output.txt"
if os.path.exists(scan_path):
    with open(scan_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            tag, ddg_str = parts[0], parts[1]
            m = re.match(r"^([A-Z]{3})([A-Z])(\d+)([A-Z])$", tag)
            if not m:
                continue
            wt3, chain, resnum, mut1 = m.groups()
            wt1 = THREE_TO_ONE.get(wt3, "?")
            res_id = int(resnum)
            ddg = float(ddg_str)
            if wt1 == mut1 or res_id in FUNCTIONAL_MASK or ddg >= -0.8:
                continue

            # Penalizza Proline/Glicine dentro le alfa-eliche o Cisteine libere (rischio ossidazione)
            biophys_warning = ""
            penalty = 0.0
            if mut1 == "C":
                biophys_warning = "[Evitare CYS libera]"
                penalty = 2.0
            elif mut1 in ("P", "G") and res_id in HELIX_CORE_RESIDS:
                biophys_warning = "[Helix-breaker in alfa-elica]"
                penalty = 2.5
            elif mut1 == "P" and res_id == 60:
                biophys_warning = "[Cautela: PRO irrigidisce loop]"
                penalty = 1.5

            qi = qi_map.get(res_id, 0.75)
            cons_label, cons_bonus = CONSENSUS_INFO.get((res_id, mut1), ("", 0.0))
            terms = energy_terms.get((res_id, mut1), {})

            # Bonus di conferma dai tool online di Alessia (se presenti i CSV in 03_mutations/)
            ddg_mpnn = mpnn_map.get((res_id, mut1))
            ddg_fp = fireprot_map.get((res_id, mut1))
            ddg_dyn = dynamut_map.get((res_id, mut1))
            cloud_bonus = 0.0
            if ddg_mpnn is not None and ddg_mpnn < 0.0:
                cloud_bonus += min(1.5, abs(ddg_mpnn))
            if ddg_fp is not None and ddg_fp < 0.0:
                cloud_bonus += min(1.5, abs(ddg_fp))
            if ddg_dyn is not None and ddg_dyn < 0.0:
                cloud_bonus += min(1.0, abs(ddg_dyn))

            # Score composito (piu' negativo = migliore candidato per MD):
            # combina ddG_FoldX + bonus vulnerabilita' MD (1 - Q_i) + bonus Consenso Ankyrin + Cloud - penalita' biofisiche
            md_vuln_bonus = max(0.0, (0.75 - qi) * 4.0)
            composite_score = ddg - md_vuln_bonus - cons_bonus - cloud_bonus + penalty

            records.append({
                "mut": f"{wt1}{res_id}{mut1}",
                "foldx_code": f"{wt1}{chain}{res_id}{mut1}",
                "res_id": res_id,
                "wt1": wt1,
                "mut1": mut1,
                "repeat": "ANK1" if res_id <= 34 else ("ANK2" if res_id <= 67 else ("ANK3" if res_id <= 100 else "ANK4")),
                "ddg_foldx": ddg,
                "qi_400K": qi,
                "rmsf_400K": rmsf_map.get(res_id, 0.0),
                "ddg_mpnn": ddg_mpnn,
                "ddg_fireprot": ddg_fp,
                "ddg_dynamut": ddg_dyn,
                "dHBond": terms.get("dHBond", 0.0),
                "dVdW_SolvH": terms.get("dVdW", 0.0) + terms.get("dSolvH", 0.0),
                "dClash": terms.get("dClash", 0.0),
                "cons_label": (cons_label + " " + biophys_warning).strip(),
                "score": composite_score
            })

records.sort(key=lambda r: r["score"])

# Aggiorna automaticamente 03_mutations/dynamut2_mutations_list.txt con le Top 30 mutazioni "pulite"
clean_candidates = [r for r in records if "[" not in r["cons_label"]]
dynamut_list_path = "03_mutations/dynamut2_mutations_list.txt"
with open(dynamut_list_path, "w") as f_dyn:
    for r in clean_candidates[:30]:
        f_dyn.write(f"A {r['mut']}\n")

out_master = "03_mutations/consensus_master_table.txt"
with open(out_master, "w") as f:
    f.write("==========================================================================================================================================================\n")
    f.write(" MATRICE DI CONSENSO MULTI-TOOL (ORTHOGONAL FUNNEL: MD 400K + FoldX Energy Breakdown + Ankyrin Consensus + ThermoMPNN/FireProt/DynaMut2)\n")
    f.write("==========================================================================================================================================================\n")
    f.write(f"{'Rank':<5} {'Mut':<8} {'Repeat':<7} {'ddG_FoldX':<10} {'<Q_i>_400K':<11} {'d(VdW+SolvH)':<13} {'dClash/Tors':<12} {'MPNN':<8} {'FireProt':<9} {'DynaMut2':<9} {'Score':<8} {'Razionale Biofisico & Consenso'}\n")
    f.write("-" * 154 + "\n")
    for idx, r in enumerate(records[:30], 1):
        mpnn_s = f"{r['ddg_mpnn']:+.2f}" if r["ddg_mpnn"] is not None else "  -  "
        fp_s = f"{r['ddg_fireprot']:+.2f}" if r["ddg_fireprot"] is not None else "   -   "
        dyn_s = f"{r['ddg_dynamut']:+.2f}" if r["ddg_dynamut"] is not None else "   -   "
        f.write(f"{idx:<5d} {r['mut']:<8s} {r['repeat']:<7s} {r['ddg_foldx']:<+10.3f} {r['qi_400K']:<11.3f} "
                f"{r['dVdW_SolvH']:<+13.3f} {r['dClash']:<+12.3f} {mpnn_s:<8s} {fp_s:<9s} {dyn_s:<9s} {r['score']:<+8.2f} {r['cons_label']}\n")

# Grafico PDF della scomposizione energetica (Energy Breakdown) per le Top 10 mutazioni "pulite" (senza CYS/PRO-breaker)
clean_top = [r for r in records if "[" not in r["cons_label"]][:10][::-1]

fig, ax = plt.subplots(figsize=(11, 5.2))
y_pos = np.arange(len(clean_top))
labels = [f"{r['mut']} ({r['repeat']})" for r in clean_top]

tot_ddg = [r["ddg_foldx"] for r in clean_top]
vdw_solv = [r["dVdW_SolvH"] for r in clean_top]
clash_tors = [r["dClash"] for r in clean_top]
hbond = [r["dHBond"] for r in clean_top]

bar_w = 0.22
ax.barh(y_pos + 1.5*bar_w, tot_ddg, height=bar_w, color="#1f77b4", edgecolor="#111111", label="ΔΔG Totale FoldX (kcal/mol)")
ax.barh(y_pos + 0.5*bar_w, vdw_solv, height=bar_w, color="#2ca02c", edgecolor="#111111", label="Δ(VdW + Solvatazione Idrofobica)")
ax.barh(y_pos - 0.5*bar_w, clash_tors, height=bar_w, color="#ff7f0e", edgecolor="#111111", label="Δ(Rilascio Clash Sterici & Torsionali)")
ax.barh(y_pos - 1.5*bar_w, hbond, height=bar_w, color="#9467bd", edgecolor="#111111", label="Δ(Legami a Idrogeno H-Bond)")

ax.axvline(0.0, color="#222222", lw=1.0)
ax.axvline(-1.0, color="#cc0000", lw=1.1, ls="--", label="Soglia -1.0 kcal/mol")
ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
ax.set_xlabel("Contributo Energetico ΔΔG (kcal/mol)  [Negativo = Stabilizzante]", fontsize=10, fontweight="bold")
ax.set_title("Scomposizione Fisico-Chimica (FoldX Energy Breakdown) dei Migliori Candidati Multi-Tool", fontsize=11.5, fontweight="bold")
ax.grid(axis="x", linestyle=":", alpha=0.6)
ax.legend(loc="lower left", frameon=True, fontsize=8.5)

plt.tight_layout()
pdf_out = "06_presentation/fig8_consensus_energy_breakdown.pdf"
plt.savefig(pdf_out, format="pdf", bbox_inches="tight")
plt.close()

print("==================================================================")
print(f" Tabella Master di Consenso salvata in : {out_master}")
print(f" Grafico Scomposizione PDF salvato in  : {pdf_out}")
print("==================================================================")
