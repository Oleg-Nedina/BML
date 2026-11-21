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
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "H1S": "H", "H2S": "H", "H3S": "H"
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items() if not k.startswith("H1") and not k.startswith("H2") and not k.startswith("H3")}

FUNCTIONAL_MASK = {
    4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36,
    41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100
}
TERMINAL_RESIDS = {1, 2, 3, 116, 117, 118}

# Residui in alfa-elica interna/esterna (HELIX record di 2MYO.pdb + Yang et al. 1998) dove PRO/GLY rompono l'elica
HELIX_CORE_RESIDS = (
    set(range(4, 11)) | set(range(16, 22)) |
    set(range(40, 46)) | set(range(50, 56)) |
    set(range(72, 78)) | set(range(83, 88)) |
    set(range(103, 111))
)

CONSENSUS_INFO = {
    (59, "G"):  ("★ Consenso ANK (GADVNA pos 25)", +2.0),
    (68, "G"):  ("★ Consenso ANK & FireProt BTC (Beta-Turn pos 1)", +2.0),
    (68, "N"):  ("★ Consenso ANK (Beta-Turn pos 1)", +1.8),
    (81, "L"):  ("★ Packing Loop Inter-Elica ANK3", +1.5),
    (103, "F"): ("★ Core Hydrophobic C-cap ANK4", +1.8),
    (103, "L"): ("★ Consenso ANK (TPLH pos 5 C-cap ANK4)", +1.8),
    (108, "L"): ("★ Packing Elica C-cap ANK4", +1.5),
    (108, "A"): ("★ Consenso ANK & FireProt BTC (Elica 4A C-cap)", +1.5),
    (105, "F"): ("★ Variante Umana Naturale (3AAA)", +1.8),
    (37, "L"):  ("★ Stabilizzazione Ingresso TPLH ANK2", +1.5),
    (37, "I"):  ("★ Stabilizzazione Ingresso TPLH ANK2", +1.2),
    (37, "F"):  ("★ Stabilizzazione Ingresso TPLH ANK2", +1.2),
    (37, "T"):  ("★ Consenso ANK (TPLH pos 3 ANK2)", +1.5),
    (37, "H"):  ("★ FireProt HighCombined (Ingresso TPLH ANK2)", +1.2),
    (60, "L"):  ("★ Rilascio Tensione Loop GADVNA ANK2", +1.2),
    (64, "V"):  ("★ Rimuove Outlier Ramachandran PRO 64", +1.5),
    (93, "G"):  ("★ FireProt HighCombined (Rilassa Loop GADKTV ANK3)", +1.5),
    (27, "F"):  ("★ Packing Loop GADVNA in ANK1", +1.0),
    (74, "A"):  ("★ Consenso ANK & FireProt BTC (TPLHLAA ANK3)", +1.5),
    (74, "I"):  ("★ Packing Core Elica 3A in ANK3", +1.0),
    (76, "I"):  ("★ FireProt BTC (Core Elica 3A in ANK3)", +1.2),
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

# 2. Carica ThermoMPNN / FireProt / DynaMut2 dai CSV di Alessia in 03_mutations/
mpnn_map = {}
if os.path.exists("03_mutations/thermompnn_2MYO.csv"):
    with open("03_mutations/thermompnn_2MYO.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ddg_val = (
                row.get("ddG (kcal/mol)") or row.get("ddG_pred")
                or row.get("ddG") or row.get("ddg")
            )
            pos_str = (row.get("Position") or row.get("position") or "").strip()
            mut_aa = (row.get("Mutant AA") or row.get("mut_aa") or "").strip()
            mut_str = (row.get("Mutation") or row.get("mutation") or row.get("mut") or "").strip()
            if pos_str.isdigit() and len(mut_aa) == 1 and ddg_val is not None:
                try:
                    mpnn_map[(int(pos_str), mut_aa.upper())] = float(ddg_val)
                except ValueError:
                    pass
            else:
                m = re.search(r"([A-Z])(\d+)([A-Z])$", mut_str)
                if m and ddg_val is not None:
                    try:
                        mpnn_map[(int(m.group(2)), m.group(3))] = float(ddg_val)
                    except ValueError:
                        pass

fireprot_map = {}
fireprot_meta = {}
if os.path.exists("03_mutations/fireprot_results.csv"):
    with open("03_mutations/fireprot_results.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pos_str = (row.get("position") or row.get("Position") or row.get("res_num") or "").strip()
            sub_str = (
                row.get("mutation") or row.get("Substitution")
                or row.get("mutant") or row.get("mut_aa") or ""
            ).strip()
            mut_str = (row.get("Mutation") or "").strip()

            res_num, mut_aa = None, None
            if pos_str.isdigit() and len(sub_str) == 1:
                res_num, mut_aa = int(pos_str), sub_str.upper()
            else:
                m = re.search(r"([A-Z])(\d+)([A-Z])", mut_str or sub_str)
                if m:
                    res_num, mut_aa = int(m.group(2)), m.group(3)

            if res_num is not None and mut_aa is not None:
                candidates_ddg = [
                    row.get("rosetta"), row.get("Rosetta"), row.get("rosetta_ddg"), row.get("ddG_Rosetta"),
                    row.get("foldx"), row.get("FoldX"), row.get("foldx_ddg"), row.get("ddG")
                ]
                ddg_num = None
                for cval in candidates_ddg:
                    if cval not in (None, "", "-", "null", "NULL", "NaN"):
                        try:
                            ddg_num = float(cval)
                            break
                        except ValueError:
                            continue

                is_btc = str(row.get("isBTC", "")).lower() == "true"
                is_he = str(row.get("isHighEnergy", "")).lower() == "true"
                is_hc = str(row.get("isHighCombined", "")).lower() == "true"
                is_maj = str(row.get("majority", "")).lower() == "true"
                is_rat = str(row.get("ratio", "")).lower() == "true"
                cons_sc = row.get("conservation_score", "")

                if ddg_num is not None:
                    fireprot_map[(res_num, mut_aa)] = ddg_num
                elif is_btc or is_he or is_hc:
                    fireprot_map[(res_num, mut_aa)] = -1.0

                fireprot_meta[(res_num, mut_aa)] = {
                    "ddg": ddg_num,
                    "isBTC": is_btc,
                    "isHighEnergy": is_he,
                    "isHighCombined": is_hc,
                    "majority": is_maj,
                    "ratio": is_rat,
                    "conservation": cons_sc
                }

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
            m = re.match(r"^([A-Z0-9]{3})([A-Z])(\d+)([A-Z])$", tag)
            if not m:
                continue
            wt3, chain, resnum, mut1 = m.groups()
            wt1 = THREE_TO_ONE.get(wt3, "?")
            res_id = int(resnum)
            ddg = float(ddg_str)
            if wt1 == mut1 or res_id in FUNCTIONAL_MASK or res_id in TERMINAL_RESIDS or ddg >= -0.8:
                continue

            # Penalizza Proline/Glicine dentro le alfa-eliche, Cisteine libere o rottura del ponte H inter-repeat HIS40->SER74
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
            elif res_id == 74 and mut1 not in ("T", "S"):
                biophys_warning = "[Cautela: SER74 OG accetta H-bond da HIS40]"
                penalty = 1.6

            qi = qi_map.get(res_id, 0.75)
            cons_label, cons_bonus = CONSENSUS_INFO.get((res_id, mut1), ("", 0.0))
            terms = energy_terms.get((res_id, mut1), {})

            # Dati dai tool online di Alessia (ThermoMPNN, FireProt 2.0, DynaMut2)
            ddg_mpnn = mpnn_map.get((res_id, mut1))
            ddg_fp = fireprot_map.get((res_id, mut1))
            fp_info = fireprot_meta.get((res_id, mut1), {})
            ddg_dyn = dynamut_map.get((res_id, mut1))

            is_btc = fp_info.get("isBTC", False)
            is_he = fp_info.get("isHighEnergy", False)
            is_hc = fp_info.get("isHighCombined", False)

            # Conteggio oggettivo dei 5 Pilastri indipendenti:
            # 1) Vulnerabilita' MD a 400K (<Qi>_400K < 0.72)
            # 2) Stabilizzazione FoldX locale (ddG < -1.0 kcal/mol)
            # 3) Fitness Strutturale AI ThermoMPNN-D (ddG_MPNN < 0.0 kcal/mol)
            # 4) FireProt 2.0 (ddG_FireProt < -0.5 oppure isBTC / isHighCombined / isHighEnergy)
            # 5) Consenso Ankyrin / Back-to-Consensus (Mosavi 2002 o FireProt isBTC)
            pillars = 0
            if qi < 0.72:
                pillars += 1
            if ddg < -1.0:
                pillars += 1
            if ddg_mpnn is not None and ddg_mpnn < 0.0:
                pillars += 1
            if (ddg_fp is not None and ddg_fp < -0.5) or is_btc or is_hc or is_he:
                pillars += 1
            if cons_bonus > 0.0 or is_btc:
                pillars += 1

            cloud_bonus = 0.0
            if ddg_mpnn is not None and ddg_mpnn < 0.0:
                cloud_bonus += min(1.5, abs(ddg_mpnn) * 2.0)
            if ddg_fp is not None and ddg_fp < 0.0:
                cloud_bonus += min(1.5, abs(ddg_fp) * 0.6)
            if is_btc:
                cloud_bonus += 1.2
            if is_hc or is_he:
                cloud_bonus += 0.5
            if ddg_dyn is not None and ddg_dyn < 0.0:
                cloud_bonus += min(1.0, abs(ddg_dyn))

            # Score composito (piu' negativo = migliore candidato per MD):
            md_vuln_bonus = max(0.0, (0.75 - qi) * 4.0)
            composite_score = ddg - md_vuln_bonus - cons_bonus - cloud_bonus + penalty

            fp_tag = ""
            if is_btc:
                fp_tag = "BTC"
            elif is_hc and is_he:
                fp_tag = "HE,HC"
            elif is_hc:
                fp_tag = "HC"
            elif is_he:
                fp_tag = "HE"

            records.append({
                "mut": f"{wt1}{res_id}{mut1}",
                "foldx_code": f"{wt1}{chain}{res_id}{mut1}",
                "res_id": res_id,
                "wt1": wt1,
                "mut1": mut1,
                "repeat": "ANK1" if res_id <= 34 else ("ANK2" if res_id <= 67 else ("ANK3" if res_id <= 100 else "ANK4")),
                "pillars": pillars,
                "ddg_foldx": ddg,
                "qi_400K": qi,
                "rmsf_400K": rmsf_map.get(res_id, 0.0),
                "ddg_mpnn": ddg_mpnn,
                "ddg_fireprot": ddg_fp,
                "fp_tag": fp_tag,
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

# Estrai anche le migliori coppie epistatiche/sinergiche da 03_mutations/thermompnn_double_2MYO.csv
double_synergies = []
if os.path.exists("03_mutations/thermompnn_double_2MYO.csv"):
    single_foldx = {(r["res_id"], r["mut1"]): r["ddg_foldx"] for r in records}
    with open("03_mutations/thermompnn_double_2MYO.csv") as f_dbl:
        for row in csv.DictReader(f_dbl):
            p1, p2 = int(row["Position 1"]), int(row["Position 2"])
            m1, m2 = row["Mutant AA 1"], row["Mutant AA 2"]
            w1, w2 = row["Wildtype AA 1"], row["Wildtype AA 2"]
            if p1 in (FUNCTIONAL_MASK | TERMINAL_RESIDS) or p2 in (FUNCTIONAL_MASK | TERMINAL_RESIDS):
                continue
            if "C" in (m1, m2):
                continue
            if (m1 in ("P", "G") and p1 in HELIX_CORE_RESIDS) or (m2 in ("P", "G") and p2 in HELIX_CORE_RESIDS):
                continue
            if (p1, m1) in single_foldx and (p2, m2) in single_foldx:
                ddg_d = float(row["ddG (kcal/mol)"])
                fx_sum = single_foldx[(p1, m1)] + single_foldx[(p2, m2)]
                if ddg_d < -0.25 and fx_sum < -2.5:
                    double_synergies.append((
                        ddg_d, fx_sum, f"{w1}{p1}{m1} + {w2}{p2}{m2}", float(row["CA-CA Distance"]),
                        mpnn_map.get((p1, m1), 0.0) + mpnn_map.get((p2, m2), 0.0)
                    ))
    double_synergies.sort(key=lambda x: x[0] + 0.25 * x[1])

out_master = "03_mutations/consensus_master_table.txt"
with open(out_master, "w") as f:
    f.write("========================================================================================================================================================================\n")
    f.write(" MATRICE DI CONSENSO MULTI-TOOL (ORTHOGONAL FUNNEL: MD 400K + FoldX 88/88 Energy Breakdown + Ankyrin Consensus + FireProt 2.0 + ThermoMPNN-D)\n")
    f.write("========================================================================================================================================================================\n")
    f.write(f"{'Rank':<5} {'Mut':<8} {'Repeat':<7} {'Pillars':<8} {'ddG_FoldX':<10} {'<Q_i>_400K':<11} {'d(VdW+SolvH)':<13} {'dClash/Tors':<12} {'MPNN':<8} {'FireProt':<13} {'DynaMut2':<9} {'Score':<8} {'Razionale Biofisico & Consenso'}\n")
    f.write("-" * 168 + "\n")
    for idx, r in enumerate(records[:35], 1):
        mpnn_s = f"{r['ddg_mpnn']:+.2f}" if r["ddg_mpnn"] is not None else "  -  "
        if r["ddg_fireprot"] is not None:
            fp_s = f"{r['ddg_fireprot']:+.2f}" + (f"({r['fp_tag']})" if r["fp_tag"] else "")
        else:
            fp_s = "   -   "
        dyn_s = f"{r['ddg_dynamut']:+.2f}" if r["ddg_dynamut"] is not None else "   -   "
        pil_s = f"{r['pillars']}/5"
        f.write(f"{idx:<5d} {r['mut']:<8s} {r['repeat']:<7s} {pil_s:<8s} {r['ddg_foldx']:<+10.3f} {r['qi_400K']:<11.4f} "
                f"{r['dVdW_SolvH']:<+13.3f} {r['dClash']:<+12.3f} {mpnn_s:<8s} {fp_s:<13s} {dyn_s:<9s} {r['score']:<+8.2f} {r['cons_label']}\n")

    if double_synergies:
        f.write("\n========================================================================================================================================================================\n")
        f.write(" TOP COPPIE DI DOPPI MUTANTI LOCALI VALUTATE DA THERMOMPNN-D (thermompnn_double_2MYO.csv) CON FOLDX CONCORDE\n")
        f.write("========================================================================================================================================================================\n")
        f.write(f"{'Coppia Doppia':<18} {'d(CA-CA)':<10} {'ddG_MPNN_Double':<17} {'Somma_MPNN_Singole':<20} {'Epistasi_MPNN':<15} {'Somma_FoldX_Singole'}\n")
        f.write("-" * 110 + "\n")
        seen_pairs = set()
        for ddg_d, fx_sum, pair_lbl, dist_ca, mp_sum in double_synergies:
            pos_pair = tuple(sorted(re.findall(r"\d+", pair_lbl)))
            if pos_pair in seen_pairs:
                continue
            seen_pairs.add(pos_pair)
            epist = ddg_d - mp_sum
            f.write(f"{pair_lbl:<18s} {dist_ca:<10.2f} {ddg_d:<+17.3f} {mp_sum:<+20.3f} {epist:<+15.3f} {fx_sum:<+15.3f}\n")
            if len(seen_pairs) >= 10:
                break

# Grafico PDF della scomposizione energetica (Energy Breakdown) per le Top 10 mutazioni "pulite" (senza CYS/PRO-breaker/H-bond breaker)
clean_top = [r for r in records if "[" not in r["cons_label"]][:10][::-1]

fig, ax = plt.subplots(figsize=(11.5, 5.5))
y_pos = np.arange(len(clean_top))
labels = [f"{r['mut']} ({r['repeat']}, {r['pillars']}/5★)" for r in clean_top]

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
ax.set_title("Scomposizione Fisico-Chimica (FoldX Energy Breakdown) dei Migliori Candidati Multi-Tool (5 Pilastri)", fontsize=11.5, fontweight="bold")
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
