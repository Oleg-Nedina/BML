#!/usr/bin/env python3
"""
07_run_foldx_buildmodel.py
Punto 3 della Scaletta del Prof. Gautieri (Laboratory_Project_2026.pdf, Slide 18-19)
& Livello 1 della Politica di Validazione (ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md):

1. Prepara la cartella dedicata 03_mutations/build_model/ con 2MYO_Repair.pdb
   e il file individual_list.txt contenente i candidati finalisti:
   - Mutanti singoli: K37L (ANK2, 5/5), T103L (ANK4, 5/5), H68G (ANK3, #1 BTC FireProt),
     D93G (ANK3, 4/5 HighC), T103F (ANK4, 4/5), A59G (ANK2), D60L (ANK2)
   - Doppi mutanti inter-repeat: K37L+T103L, K37L+H68G, H68G+T103L, K37L+D93G
   - Doppio mutante locale sinergico (ThermoMPNN-D): A59G+D60L
   - Triplo mutante multi-repeat (ANK2+ANK3+ANK4): K37L+H68G+T103L
2. Esegue FoldX 5 BuildModel con numberOfRuns=5 (generando sia i mutanti 2MYO_Repair_<i>_<run>.pdb
   sia i "matched WT references" WT_2MYO_Repair_<i>_<run>.pdb richiesti dalla Slide 19!).
3. Analizza Average_2MYO_Repair.fxout e Dif_2MYO_Repair.fxout calcolando:
   - Media ± Deviazione Standard (μ ± σ) su 5 run di ΔΔG_BuildModel e dei singoli termini fisici
   - Ciclo energetico di Additività vs Epistasi:
     ΔΔΔG_coupling = ΔΔG_Combined - Σ ΔΔG_Singles
4. Seleziona e copia nella cartella pulita 03_mutations/final_models_pdb/ il conformero a più
   bassa energia ( e il relativo matched WT ) pronto per la solvatazione del Punto 4.
5. Esporta il report 03_mutations/buildmodel_epistasis_report.txt e il grafico vettoriale
   06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf.
"""

import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt

BM_DIR = "03_mutations/build_model"
FINAL_PDB_DIR = "03_mutations/final_models_pdb"
os.makedirs(BM_DIR, exist_ok=True)
os.makedirs(FINAL_PDB_DIR, exist_ok=True)
os.makedirs("06_presentation", exist_ok=True)

# 1. Lista ordinata dei mutanti da costruire con FoldX BuildModel (numberOfRuns=5)
# Formato: (Etichetta_Umana, Codice_individual_list, Categoria, Componenti_Singole, Descrizione)
VARIANTS = [
    ("K37L",            "KA37L;",               "Single", ["K37L"],                 "ANK2 (Ingresso TPLH, 5/5 Pilastri, <Q_37>=0.485)"),
    ("T103L",           "TA103L;",              "Single", ["T103L"],                "ANK4 (C-cap Core Consenso TPLH, 5/5 Pilastri, <Q_103>=0.508)"),
    ("H68G",            "HA68G;",               "Single", ["H68G"],                 "ANK3 (Beta-Turn pos 1, #1 FireProt BTC, <Q_68>=0.714)"),
    ("D93G",            "DA93G;",               "Single", ["D93G"],                 "ANK3 (Loop GADKTV, 4/5 Pilastri FireProt HighC, <Q_93>=0.543)"),
    ("T103F",           "TA103F;",              "Single", ["T103F"],                "ANK4 (C-cap Core Aromatico, 4/5 Pilastri, <Q_103>=0.508)"),
    ("A59G",            "AA59G;",               "Single", ["A59G"],                 "ANK2 (Ingresso Loop GADVNA Consenso, <Q_59>=0.640)"),
    ("D60L",            "DA60L;",               "Single", ["D60L"],                 "ANK2 (Loop GADVNA Rilascio Tensione, <Q_60>=0.465)"),
    ("K37L+T103L",      "KA37L,TA103L;",        "Double", ["K37L", "T103L"],        "Doppio A: ANK2 + ANK4 (5/5 + 5/5 Pilastri, Poli Opposti)"),
    ("K37L+H68G",       "KA37L,HA68G;",         "Double", ["K37L", "H68G"],         "Doppio B: ANK2 + ANK3 (5/5 + #1 FireProt BTC)"),
    ("H68G+T103L",      "HA68G,TA103L;",        "Double", ["H68G", "T103L"],        "Doppio C: ANK3 + ANK4 (#1 FireProt BTC + 5/5 Pilastri)"),
    ("K37L+D93G",       "KA37L,DA93G;",         "Double", ["K37L", "D93G"],         "Doppio D: ANK2 + ANK3 (5/5 + FireProt HighCombined)"),
    ("A59G+D60L",       "AA59G,DA60L;",         "Double", ["A59G", "D60L"],         "Doppio E: ANK2 Locale Sinergico (ThermoMPNN-D Epistatico)"),
    ("K37L+H68G+T103L", "KA37L,HA68G,TA103L;",  "Triple", ["K37L", "H68G", "T103L"],"Triplo Consenso: ANK2 + ANK3 + ANK4 (1 Bullone per Repeat)"),
]

N_RUNS = 5

# Copia 2MYO_Repair.pdb nella cartella build_model
shutil.copyfile("03_mutations/2MYO_Repair.pdb", os.path.join(BM_DIR, "2MYO_Repair.pdb"))

# Scrivi individual_list.txt (sia in build_model/ sia in 03_mutations/ per riferimento)
ind_list_bm = os.path.join(BM_DIR, "individual_list.txt")
ind_list_root = "03_mutations/individual_list.txt"
with open(ind_list_bm, "w") as f_bm, open(ind_list_root, "w") as f_rt:
    for _, code, _, _, _ in VARIANTS:
        f_bm.write(code + "\n")
        f_rt.write(code + "\n")

# Scrivi build_model.cfg
cfg_bm = os.path.join(BM_DIR, "build_model.cfg")
with open(cfg_bm, "w") as f_cfg:
    f_cfg.write("command=BuildModel\n")
    f_cfg.write("pdb=2MYO_Repair.pdb\n")
    f_cfg.write("mutant-file=individual_list.txt\n")
    f_cfg.write(f"numberOfRuns={N_RUNS}\n")

print("==========================================================================")
print(f" AVVIO FOLDX 5 BUILDMODEL ({len(VARIANTS)} Varianti x {N_RUNS} Run = {len(VARIANTS)*N_RUNS} Modelli + Matched WT)")
print("==========================================================================")

# Pulisci eventuali vecchi file .fxout in BM_DIR prima di lanciare
for old_fx in ("Average_2MYO_Repair.fxout", "Dif_2MYO_Repair.fxout", "Raw_2MYO_Repair.fxout", "PdbList_2MYO_Repair.fxout"):
    p = os.path.join(BM_DIR, old_fx)
    if os.path.exists(p):
        os.remove(p)

with open(os.path.join(BM_DIR, "build_model.log"), "w") as flog:
    subprocess.run(["foldx", "-f", "build_model.cfg"], cwd=BM_DIR, stdout=flog, stderr=subprocess.STDOUT, check=True)

print(" -> FoldX BuildModel terminato con successo! Analisi di Dif_2MYO_Repair.fxout...")

# 2. Parsing di Dif_2MYO_Repair.fxout (contiene le 5 run per ciascuna delle 13 varianti = 65 righe)
dif_path = os.path.join(BM_DIR, "Dif_2MYO_Repair.fxout")
dif_rows = []
with open(dif_path) as f:
    for line in f:
        line = line.strip()
        if not line or not line.startswith("2MYO_Repair_"):
            continue
        parts = line.split("\t")
        if len(parts) < 16:
            continue
        pdb_name = parts[0]
        vals = [float(x) for x in parts[1:]]
        # Indici colonne in Dif_2MYO_Repair.fxout:
        # 0:total, 1:BackHbond, 2:SideHbond, 3:VdW, 4:Elec, 5:SolvP, 6:SolvH,
        # 7:VdWClash, 8:EntropySide, 9:EntropyMain, 13:TorsionalClash, 14:BackboneClash
        dif_rows.append({
            "pdb": pdb_name,
            "ddG": vals[0],
            "dHBond": vals[1] + vals[2],
            "dVdW": vals[3],
            "dElec": vals[4],
            "dSolvP": vals[5],
            "dSolvH": vals[6],
            "dVdW_SolvH": vals[3] + vals[6],
            "dClash": vals[7] + vals[13] + vals[14],
            "dEntropy": vals[8] + vals[9],
        })

assert len(dif_rows) == len(VARIANTS) * N_RUNS, (
    f"Attese {len(VARIANTS)*N_RUNS} righe in Dif_2MYO_Repair.fxout, trovate {len(dif_rows)}"
)

# Copia i file .fxout principali anche in 03_mutations/ per comodità
for fx_name in ("Average_2MYO_Repair.fxout", "Dif_2MYO_Repair.fxout", "Raw_2MYO_Repair.fxout"):
    shutil.copyfile(os.path.join(BM_DIR, fx_name), os.path.join("03_mutations", fx_name))

# 3. Raggruppa le N_RUNS per ciascuna variante e seleziona il conformero migliore (min ddG)
results = {}
for idx, (label, code, cat, comps, desc) in enumerate(VARIANTS, 1):
    chunk = dif_rows[(idx - 1) * N_RUNS : idx * N_RUNS]
    ddgs = [r["ddG"] for r in chunk]
    hbonds = [r["dHBond"] for r in chunk]
    vdw_solvh = [r["dVdW_SolvH"] for r in chunk]
    clashes = [r["dClash"] for r in chunk]
    entropies = [r["dEntropy"] for r in chunk]

    best_run_idx = int(np.argmin(ddgs))  # 0..4
    best_pdb_name = f"2MYO_Repair_{idx}_{best_run_idx + 1}.pdb"
    best_wt_name = f"WT_2MYO_Repair_{idx}_{best_run_idx + 1}.pdb"

    # Copia il conformero migliore e il suo matched WT in 03_mutations/final_models_pdb/
    safe_lbl = label.replace("+", "_")
    src_mut = os.path.join(BM_DIR, best_pdb_name)
    src_wt = os.path.join(BM_DIR, best_wt_name)
    if os.path.exists(src_mut):
        shutil.copyfile(src_mut, os.path.join(FINAL_PDB_DIR, f"2myo_mut_{safe_lbl}.pdb"))
    if os.path.exists(src_wt):
        shutil.copyfile(src_wt, os.path.join(FINAL_PDB_DIR, f"2myo_matched_WT_for_{safe_lbl}.pdb"))

    results[label] = {
        "idx": idx,
        "label": label,
        "code": code,
        "cat": cat,
        "comps": comps,
        "desc": desc,
        "ddg_mean": float(np.mean(ddgs)),
        "ddg_sd": float(np.std(ddgs, ddof=1)) if N_RUNS > 1 else 0.0,
        "ddg_min": float(np.min(ddgs)),
        "dHBond_mean": float(np.mean(hbonds)),
        "dVdW_SolvH_mean": float(np.mean(vdw_solvh)),
        "dClash_mean": float(np.mean(clashes)),
        "dEntropy_mean": float(np.mean(entropies)),
        "best_pdb": f"03_mutations/final_models_pdb/2myo_mut_{safe_lbl}.pdb",
        "best_wt": f"03_mutations/final_models_pdb/2myo_matched_WT_for_{safe_lbl}.pdb",
    }

# 4. Calcolo del Ciclo di Additività vs Epistasi Statica per i Doppi/Tripli Mutanti
for label, info in results.items():
    if info["cat"] in ("Double", "Triple"):
        sum_singles = sum(results[c]["ddg_mean"] for c in info["comps"])
        coupling = info["ddg_mean"] - sum_singles
        info["sum_singles"] = sum_singles
        info["coupling"] = coupling
    else:
        info["sum_singles"] = info["ddg_mean"]
        info["coupling"] = 0.0

# 5. Scrittura del Report Completo 03_mutations/buildmodel_epistasis_report.txt
report_path = "03_mutations/buildmodel_epistasis_report.txt"
with open(report_path, "w") as f:
    f.write("========================================================================================================================================================\n")
    f.write(" REPORT LIVELLO 1 (PUNTO 3 SCALETTA PROF — SLIDE 19): FOLDX 5 BUILDMODEL (n=5 RUNS CON MATCHED WT REFERENCES) & ANALISI DI EPISTASI\n")
    f.write("========================================================================================================================================================\n\n")
    f.write("1. RISULTATI ENERGETICI BUILDMODEL (MEDIA ± SD SU 5 RUN INDIPENDENTI CONTRO MATCHED WT_2MYO_Repair_*.pdb)\n")
    f.write("-" * 152 + "\n")
    f.write(f"{'ID':<4} {'Variante':<17} {'Tipo':<8} {'ddG_Mean ± SD':<18} {'ddG_Best':<10} {'d(VdW+SolvH)':<14} {'dHBond':<10} {'dClash':<10} {'dEntropy':<10} {'Razionale & Modulo'}\n")
    f.write("-" * 152 + "\n")
    for label, _, _, _, _ in VARIANTS:
        r = results[label]
        ddg_s = f"{r['ddg_mean']:+.3f} ± {r['ddg_sd']:.3f}"
        f.write(f"{r['idx']:<4d} {r['label']:<17s} {r['cat']:<8s} {ddg_s:<18s} {r['ddg_min']:<+10.3f} "
                f"{r['dVdW_SolvH_mean']:<+14.3f} {r['dHBond_mean']:<+10.3f} {r['dClash_mean']:<+10.3f} {r['dEntropy_mean']:<+10.3f} {r['desc']}\n")

    f.write("\n========================================================================================================================================================\n")
    f.write("2. CICLO DI ADDITIVITA' VS EPISTASI STATICA DEI MUTANTI COMBINATI (ΔΔΔG_coupling = ΔΔG_Combined - Σ ΔΔG_Singles)\n")
    f.write("========================================================================================================================================================\n")
    f.write(f"{'Variante Combinata':<19} {'Componenti':<22} {'ddG_BuildModel (μ±σ)':<22} {'Somma_Singole (ΣddG)':<22} {'ΔΔΔG_coupling':<15} {'Verdetto Termodinamico (Gate 1)'}\n")
    f.write("-" * 152 + "\n")
    for label, _, cat, _, _ in VARIANTS:
        if cat not in ("Double", "Triple"):
            continue
        r = results[label]
        ddg_s = f"{r['ddg_mean']:+.3f} ± {r['ddg_sd']:.3f}"
        if abs(r["coupling"]) <= 0.15:
            verdict = "✅ PERFETTAMENTE ADDITIVA (Zero interferenza sterica -> Gate 1 Superato!)"
        elif r["coupling"] < -0.15:
            verdict = "🌟 SINERGIA EPISTATICA POSITIVA (Più stabile della somma delle singole!)"
        else:
            verdict = "⚠️ EPISTASI SUB-ADDITIVA (Lieve sovrapposizione locale)"
        f.write(f"{r['label']:<19s} {' + '.join(r['comps']):<22s} {ddg_s:<22s} {r['sum_singles']:<+22.3f} {r['coupling']:<+15.3f} {verdict}\n")

    f.write("\n========================================================================================================================================================\n")
    f.write("3. MODELLI 3D (.PDB) SALVATI E PRONTI PER LA SOLVATAZIONE DEL PUNTO 4 (MATCHED MD COMPARISON)\n")
    f.write("========================================================================================================================================================\n")
    for label in ("K37L", "T103L", "H68G", "K37L+T103L", "K37L+H68G", "K37L+H68G+T103L"):
        r = results[label]
        f.write(f"  • {label:<17s} -> Mutante: {r['best_pdb']:<55s} | Matched WT: {r['best_wt']}\n")

# 6. Generazione Grafico Vettoriale PDF (06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.6))
plt.subplots_adjust(wspace=0.28)

# Pannello A: ΔΔG BuildModel (μ ± σ su 5 run) per tutte le varianti
plot_order = [v[0] for v in VARIANTS][::-1]
y_pos = np.arange(len(plot_order))
means = [results[k]["ddg_mean"] for k in plot_order]
sds = [results[k]["ddg_sd"] for k in plot_order]
colors = []
for k in plot_order:
    c = results[k]["cat"]
    if c == "Single":
        colors.append("#1f77b4")
    elif c == "Double":
        colors.append("#2ca02c")
    else:
        colors.append("#9467bd")

bars = ax1.barh(y_pos, means, xerr=sds, color=colors, edgecolor="#111111", height=0.70,
                capsize=3, error_kw={"lw": 1.1, "ecolor": "#111111"})
ax1.axvline(0.0, color="#222222", lw=1.0)
ax1.axvline(-1.0, color="#cc0000", lw=1.1, ls="--", label="Soglia -1.0 kcal/mol")
for bar, m_val, sd_val in zip(bars, means, sds):
    ax1.text(m_val - sd_val - 0.08, bar.get_y() + bar.get_height() / 2.0,
             f"{m_val:+.2f}±{sd_val:.2f}", va="center", ha="right", fontsize=8, fontweight="bold")

ax1.set_yticks(y_pos)
ax1.set_yticklabels(plot_order, fontsize=9, fontweight="bold")
ax1.set_xlim(min(means) - 1.3, 0.4)
ax1.set_xlabel("FoldX BuildModel ΔΔG = ΔG_mut − ΔG_WT_matched (kcal/mol, n=5)", fontsize=9.5, fontweight="bold")
ax1.set_title("A. Stabilità FoldX BuildModel (Singole, Doppie e Tripla, μ ± σ)", fontsize=11, fontweight="bold")
ax1.grid(axis="x", linestyle=":", alpha=0.6)

# Pannello B: Confronto ΔΔG_BuildModel Reale vs Somma Teorica delle Singole (Test di Additività/Epistasi)
comb_keys = [v[0] for v in VARIANTS if v[2] in ("Double", "Triple")]
x_idx = np.arange(len(comb_keys))
w = 0.36
actual_ddg = [results[k]["ddg_mean"] for k in comb_keys]
actual_sd = [results[k]["ddg_sd"] for k in comb_keys]
sum_ddg = [results[k]["sum_singles"] for k in comb_keys]

ax2.bar(x_idx - w/2, sum_ddg, width=w, color="#a6cee3", edgecolor="#111111", label="Somma Teorica Singole (Σ ΔΔG_i)")
ax2.bar(x_idx + w/2, actual_ddg, yerr=actual_sd, width=w, color="#2ca02c", edgecolor="#111111",
        capsize=3, label="BuildModel Combinato Reale (ΔΔG_1+2)")
ax2.axhline(0.0, color="#222222", lw=1.0)
for i, k in enumerate(comb_keys):
    coup = results[k]["coupling"]
    y_txt = min(sum_ddg[i], actual_ddg[i]) - actual_sd[i] - 0.32
    ax2.text(i, y_txt, f"ΔΔΔG =\n{coup:+.2f}", ha="center", va="top", fontsize=7.8,
             fontweight="bold", color="#0b5394" if coup <= 0.15 else "#b45f06")

ax2.set_xticks(x_idx)
ax2.set_xticklabels(comb_keys, rotation=22, ha="right", fontsize=8.5, fontweight="bold")
ax2.set_ylim(min(min(actual_ddg), min(sum_ddg)) - 1.2, 0.4)
ax2.set_ylabel("ΔΔG (kcal/mol)  [Più negativo = Più stabile]", fontsize=9.5, fontweight="bold")
ax2.set_title("B. Ciclo di Additività vs Epistasi dei Mutanti Combinati", fontsize=11, fontweight="bold")
ax2.grid(axis="y", linestyle=":", alpha=0.6)
ax2.legend(loc="upper right", frameon=True, fontsize=8.5)

plt.tight_layout()
pdf_out = "06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf"
plt.savefig(pdf_out, format="pdf", bbox_inches="tight")
plt.close()

print("==========================================================================")
print(f" Report BuildModel & Epistasi salvato in : {report_path}")
print(f" Modelli 3D PDB finali salvati in        : {FINAL_PDB_DIR}/")
print(f" Grafico PDF vettoriale salvato in       : {pdf_out}")
print("==========================================================================")
