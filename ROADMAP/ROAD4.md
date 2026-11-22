# 🗺️ ROAD4.md — Fase 4: Simulazioni MD Comparative Multi-Seed ("Matched Comparison", `04_md_simulations/`)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Sistemi in Confronto Controllato (`400 K`, $n = 3$ Velocity Seeds):**
- **`WT`** (`2MYO` Wild-Type: Seed 1 già completata in `02_wt_baseline/wt_400K.dcd`; Seed 2 e 3 in `04_md_simulations/WT/`)
- **`Mut1 = K37L`** (`04_md_simulations/K37L/`, `ANK2`, 5/5 Pilastri, $\Delta\Delta G_{\text{BuildModel}} = -2.44 \pm 0.87\text{ kcal/mol}$)
- **`Mut2 = T103L`** (`04_md_simulations/T103L/`, `ANK4` C-cap, 5/5 Pilastri, $\Delta\Delta G_{\text{BuildModel}} = -1.15 \pm 0.14\text{ kcal/mol}$)
- **`Double = K37L_T103L`** (`04_md_simulations/K37L_T103L/`, `ANK2 + ANK4`, 100% Additiva $\Delta\Delta G_{\text{BuildModel}} = -3.65 \pm 0.66\text{ kcal/mol}$, `Best = -4.77 kcal/mol`)
- *(Opzionali già pronti: **`H68G`** #1 BTC FireProt in `ANK3` e **`K37L_H68G_T103L`** Triplo Mutante `-5.03 ± 0.21 kcal/mol`)*

**Stato della Fase 4:**
- ✅ **Passo 4.1 (Setup Topologico Matched `CHARMM36m + TIP3P + 0.15 M NaCl` & Gate 2 QC):** **COMPLETATO AL 100%** ([`04_md_simulations/matched_systems_qc_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/matched_systems_qc_report.txt) e [`04_md_simulations/RISULTATI_E_COMMENTO_FASE4_SETUP.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/RISULTATI_E_COMMENTO_FASE4_SETUP.md)).
- ✅ **Passo 4.2a (Screening MD `Seed 1` a `400 K` su `WT`, `K37L`, `T103L`, `K37L_T103L`, `K37L_H68G_T103L` + Plot Singoli Separati):** **COMPLETATO AL 100% (22 Nov 2026)** tramite [`04_md_simulations/run_4_screening_400K.sh`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/run_4_screening_400K.sh). Tutti i plot singoli sono in [`05_analysis/single_runs/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs) e analizzati in **[`05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md)**.
- 🔄 **Passo 4.2b (Campagna di Validazione Finale: Repliche `Seed 2 & 3` + Test con Pre-Equilibrazione Restrained a `300 K`):** **PRONTO AL LANCIO**.

---

## 📚 1. Fonti e Documentazione di Riferimento per la Fase 4

1. **Linee Guida Ufficiali del Progetto:**  
   [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 3** (*Learning Objective 3: Controlled computational comparison — Use consistent modelling and MD settings for WT and variants*), **Slide 10** (*Multiple velocity seeds*), **Slide 15** (*Step 4: Matched MD comparison*), **Slide 21** (*Computational budget and pilot benchmark*).
2. **Slide delle Lezioni del Corso:**  
   [`Lecture03 - Molecular Dynamics Setup.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture03%20-%20Molecular%20Dynamics%20Setup.pdf) (Slide 32–50: *Explicit water TIP3P, PBC, PME neutrality & counter-ions*); [`Lecture04 - Running Molecular Dynamics.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture04%20-%20Running%20Molecular%20Dynamics.pdf) (Slide 12, 15, 18–21, 26–27: *Maxwell-Boltzmann velocity seeds, NPT ensemble & reproducibility*).
3. **Paper di Riferimento:**  
   *Bekker GJ, Ma B, Kamiya N (2019).* **Protein Science**, 28(2):429–438 (pp. 433–436, Table II).

---

## ✅ 2. Riepilogo del Passo 4.1 Già Completato (Controllo di Qualità Topologico — Gate 2)

Con il comando `vmd -dispdev text -e scripts/08_build_matched_mutants_charmm36m.tcl`, avete già generato tutti i domini solvatati e verificato il superamento del **Gate 2** ([`04_md_simulations/matched_systems_qc_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/matched_systems_qc_report.txt)):
- Tutte le varianti condividono la stessa identica scatola periodica **`72.589 × 58.517 × 58.883 Å`** (`~23.070` atomi, `~7.070` molecole `TIP3P`).
- Tutti i sistemi hanno carica elettrica netta esattamente pari a **`+0.0000 e`** (`25 Na+ / 20 Cl-` per le varianti con carica proteica `-5.00 e`; `26 Na+ / 20 Cl-` per le varianti contenenti `K37L` con carica proteica `-6.00 e`).
- Per ciascuna variante sono stati creati i **3 file `.conf` indipendenti a `400 K`** (`1.000` step minimizzazione + `1.250.000` step = `2.5 ns` NPT):
  - `*_400K_seed1.conf` (`seed 12345`)
  - `*_400K_seed2.conf` (`seed 23456`)
  - `*_400K_seed3.conf` (`seed 34567`)

---

## 🚀 3. Passo 4.2 — Come Lanciare le Simulazioni MD della Fase 4 (Strategia "Seed 1 First" + Multi-Seed)

Abbiamo preparato lo script batch intelligente **[`04_md_simulations/run_batch_phase4.sh`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/run_batch_phase4.sh)** che:
- **Salta automaticamente** qualsiasi replica il cui file `.dcd` sia già stato completato (così puoi fermarlo e riprenderlo in qualsiasi momento senza mai perdere il lavoro fatto!);
- Accetta come parametro la modalità di lancio (`seed1`, `oleg`, `alessia`, `all` o `bonus`) e il numero di core fisici (`8` di default);
- Al termine di ogni batch, lancia in automatico l'analisi VMD ([`scripts/09_analyze_matched_validation.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/09_analyze_matched_validation.tcl)) e aggiorna la **Table II di Bekker et al.** e i grafici PDF ([`scripts/10_plot_final_validation_bekker.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/10_plot_final_validation_bekker.py))!

### 💻 Comandi Pronti per Lanciare i Batch (dalla cartella `BML/`):

#### 🌙 Batch 1 (Consigliato per la Prima Notte — `seed1` dei 3 mutanti principali `K37L`, `T103L`, `K37L_T103L`):
Poiché `WT Seed 1` (`02_wt_baseline/wt_400K.dcd`) è già pronta, questo batch esegue solo le 3 run `seed1` di `K37L`, `T103L` e `K37L_T103L` e a fine run compila già la prima **Table II completa su tutti e 4 i sistemi a parità di `seed 12345`**!
```bash
nohup bash 04_md_simulations/run_batch_phase4.sh seed1 8 > 04_md_simulations/phase4_seed1.log 2>&1 &
```
*(Per controllare a che punto è in tempo reale: `tail -f 04_md_simulations/phase4_seed1.log`).*

#### 🤝 Batch 2 (Divisione dei restanti `Seed 2` e `Seed 3` tra Oleg e Alessia):
- **Sul PC Linux di Oleg** (esegue `Seed 2` e `Seed 3` di `K37L` e `K37L_T103L`):
  ```bash
  nohup bash 04_md_simulations/run_batch_phase4.sh oleg 8 > 04_md_simulations/phase4_oleg.log 2>&1 &
  ```
- **Sul PC Windows di Alessia** (esegue `Seed 2` e `Seed 3` di `WT` e `T103L` — vedi guida [`ROAD_PER_ALE/03_GUIDA_FASE4_E_RISULTATI_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROAD_PER_ALE/03_GUIDA_FASE4_E_RISULTATI_CONSENSO.md)):
  ```bash
  bash 04_md_simulations/run_batch_phase4.sh alessia 8
  ```
*(Oppure, se vuoi far girare tutte le repliche rimanenti in sequenza sul tuo PC Linux: `nohup bash 04_md_simulations/run_batch_phase4.sh all 8 > 04_md_simulations/phase4_all.log 2>&1 &`).*
