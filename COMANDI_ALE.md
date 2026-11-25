# 📘 COMANDI_ALE.md — Comandi Operativi Wild-Type & Scaletta Scheletrica per la Presentazione Finale

**Corso:** Biomolecular Modelling Laboratory (2026) — Prof. Alfonso Gautieri  
**Team #10:** Alessia Rigoni & Oleg Nedina  
**Proteina:** Myotrophin (*Rattus norvegicus*, PDB: `2MYO`, 118 aa)  
**Presentazione Finale:** 📅 **Mercoledì 16 Dicembre — Max 10 minuti** ([`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf), Slide 4)

> **Nota sul workflow di coppia:** Oleg eseguirà le simulazioni multi-seed delle mutazioni e ti inoltrerà man mano i pacchetti di risultati numerici, i grafici e le considerazioni già strutturate; tu gestisci il **triplicato del Wild-Type (`WT`)**, la validazione del `WT` e l'espansione/integrazione dei contenuti nei report (`report/`) e nella presentazione finale.

---

## ⚙️ PARTE A — Comandi Operativi per Alessia (Wild-Type `WT` a `400 K`)

Esegui tutti i comandi dalla cartella principale `BML/`.

### 1. Lancio del Triplicato Wild-Type (`3 velocity seeds`: `12345`, `23456`, `34567` a `400 K`, `2.5 ns`)
Stesso identico protocollo e stessi 3 semi delle mutazioni di Oleg (resume automatico se un seed è già finito):
```bash
bash 04_md_simulations/WT/run_WT_3seeds.sh 8
```
*(Sostituisci `8` con il numero di core fisici del tuo PC; su Windows puoi usare Git Bash / WSL oppure lanciare NAMD 3 sui singoli file `04_md_simulations/WT/wt_400K_seed1.conf`, `seed2.conf`, `seed3.conf`).*

### 2. Controllo Qualità & Validazione delle Run del Wild-Type
- **Verifica completamento senza errori (`End of program` e `0` `FATAL ERROR`):**
  ```bash
  grep -E "Benchmark|End of program|FATAL ERROR" 04_md_simulations/WT/wt_400K_seed*.log
  ```
- **Verifica preparazione topologica WT (già pronta):**
  ```bash
  vmd -dispdev text -e report/02_scripts/01_validate_wt_preparation.tcl
  ```
- **Analisi completa WT (`Q-value` con `Frame 0 = 1° frame DCD` e confini *Yang et al. 1998*, `RMSD`, `RMSF`, `SASA`, `Rgyr`, `H-Bonds 3.5 Å / 30°`, `7 α-eliche H1–H7`):**
  ```bash
  vmd -dispdev text -e scripts/03_analyze_wt_baseline.tcl
  vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl
  vmd -dispdev text -e scripts/12_analyze_hbonds_and_helices.tcl
  ```
- **Rigenerazione di tutti i plot del Wild-Type (sia i tuoi in `report/02_utilities/` sia i comparativi in `06_presentation/`):**
  ```bash
  python3 report/02_scripts/02_plot_rmsd_q_wt.py
  python3 report/02_scripts/03_plot_wt_remaining.py
  python3 report/02_scripts/04_plot_wt_structural_regions.py
  python3 scripts/plot_wt_baseline.py
  python3 scripts/10_plot_final_validation_bekker.py
  python3 scripts/12_plot_hbonds_and_helices.py
  ```

### 3. Sincronizzazione Git con Oleg
```bash
git pull --no-rebase --no-edit && git push
```

---

## 🗂️ PARTE B — Scaletta Scheletrica per la Presentazione Finale (Allineata a `Laboratory_Project_2026.pdf`, Slide 3, 4, 15, 19, 21–23)

*(Struttura scheletrica: punti chiave da toccare + dove reperire i dati + script/grafici/VMD da usare).*

---

### 1️⃣ Introduzione e Analisi della Proteina Assegnata (`2MYO`)
*Requisiti `Laboratory_Project_2026.pdf`: Slide 3 (LO1: connect sequence & structure to biological function), Slide 4 (Introduction on target protein), Slide 7 (improve stability without affecting function), Slide 14–15 (chain, biological assembly, missing residues/atoms).*

- **Punti chiave da inserire:**
  - Identità, organismo (`Rattus norvegicus`), peso (~13 kDa), lunghezza (118 aa), architettura **Ankyrin-repeat** (`7 α-eliche H1–H7`, 46.6% STRIDE, 2 loop $\beta$-hairpin sporgenti `32–36` e `65–69`).
  - **Funzione biologica & basi delle ipotesi:** legame e inibizione sterica dell'**Actin-Capping Protein (CP)** (*Takeda et al., PLoS Biol. 2010*; *Yang et al., Structure 1998*) + ruolo intracellulare (*Gupta et al. 2002*; *Knuefermann et al. 2006*).
  - **Validazione sperimentale `2MYO`:** Solution NMR (*minimized average*, 297 K, pH 6.2), Monomero (`Chain A`), `118/118` residui continui (`902` atomi pesanti, `0` gap/mancanti), `3` cisteine tutte ridotte ($d_{S_\gamma\text{--}S_\gamma} \ge 9.43\text{ Å}$ $\rightarrow$ `0` ponti disolfuro).
  - **Dettagli di ricerca avanzata a supporto delle ipotesi:**
    - SASA/RSA statica (*Tien et al. 2013*): `19` residui sepolti (`RSA < 0.05`, core `ALA42, ILE50, LEU51, LEU54`), `21` H-bonds, `2` ponti salini statici (`ASP27–ARG30`, `ASP65–LYS97`).
    - `REMARK 500`: `35` outlier Ramachandran (inclusi `HIS 68`, `ASP 93`, `THR 103`) e close contact `HIS80 O – SER82 N` (`2.14 Å`) $\rightarrow$ limiti del modello NMR singolo.
    - **Functional Mask oggettiva (`3AAA`):** sovrapposizione `2MYO` su `3AAA` ($\text{RMSD} = 2.365\text{ Å}$, identità 99.2%) $\rightarrow$ **24 residui a $< 4.5\text{ Å}$ da CP resi intoccabili** (garantisce il vincolo funzionale della Slide 7).
- **Dove reperire dati e fonti:**
  - [`report/01_protein_characterization.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/01_protein_characterization.md)
  - [`SCALETTA_PROF/PUNTO_1_ANALISI_PROTEINA.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/SCALETTA_PROF/PUNTO_1_ANALISI_PROTEINA.md) & [`01_structures/RISULTATI_E_COMMENTO_FASE1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/RISULTATI_E_COMMENTO_FASE1.md)
  - Tabelle: [`report/01_utilities/02_ankyrin_region_summary.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/01_utilities/02_ankyrin_region_summary.csv), [`report/01_utilities/07_sasa_2myo_rsa.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/01_utilities/07_sasa_2myo_rsa.csv)
- **Grafici già pronti & Comandi VMD 3D:**
  - Immagini: `report/01_utilities/01_overall_protein_fold.png`, `02_ankyrin_modules.png`, `03_helix_turn_helix.png`, `04_hairpin_loops.png`, `08_rsa_distribution.png`, `09_hydrophobic_packing.png`.
  - VMD 3D (struttura + Functional Mask `3AAA`): `vmd -e scripts/viz_01_structure_and_mask.tcl`

---

### 2️⃣ Characterization and Baseline (`WT` a `300 K` vs `400 K`)
*Requisiti `Laboratory_Project_2026.pdf`: Slide 9–13 (Bekker et al. 2019, scelta 400 K e Q-value), Slide 15 (prepare WT, pilot simulation, candidate regions), Slide 21 (benchmark first & computational budget), Slide 22–23 (pbc wrap, align, RMSF `analysisMD.tcl`, Q-value `VMDextensions.tcl`).*

- **Punti chiave da inserire:**
  - **Setup sistema `CHARMM36m`:** assegnazione istidine a pH neutro su base geometrica H-bond entro $3.5\text{ Å}$ (`HIS 40, 68, 80 -> HSD`; `HIS 67 -> HSE` tra `LYS 37` e `GLY 35 O`), rimozione idrogeni NMR (`protein and noh`), `psfgen` ($\Delta C_\alpha = 0.0000\text{ Å}$), solvatazione `TIP3P` (`12 Å`, `6.535` acque), ionizzazione `0.15 M NaCl` (`23 Na+ / 18 Cl-`, `21.461` atomi, $q_{\text{tot}} = 0.0000\text{ e}$).
  - **Benchmark & Budget Computazionale (Slide 21):** `namd3 +p8` = `12.68 ns/day` (`+40%` rispetto a `+p16` SMT); calcolo esplicito di `Total MD time = systems × T × replicas × length`.
  - **Perché `400 K` (Slide 10–12):** `300 K` (`1.5 ns`) è il controllo stabile ($\text{RMSD} \approx 2.6\text{ Å}$, $R_g = 14.5\text{ Å}$); `400 K` (`2.5 ns`) accelera la perdita di contatti senza denaturazione immediata, massimizzando la correlazione con $T_m$ (*Bekker et al. 2019*: $r=0.79$ `all-all`, $r=0.84$ `hydrophilic-all`).
  - **Nota metodologica (confermata dal Prof.):** citare e giustificare la scelta di iniziare NPT dopo `1.000` step di minimizzazione (senza pre-equilibrazione restrained separata) valutando le medie sull'**ultimo 30% equilibrato (`1.75–2.50 ns`)** e usando il **`Frame 0` del `.dcd`** come riferimento nativo ($Q(0)=1.00$).
  - **Perché abbiamo diviso la proteina in 4 regioni (`Yang et al., 1998`):** `N-terminal (1–27: H1–H2)`, `Ankyrin repeat 1 (28–60: loop 32–36, H3–H4)`, `Ankyrin repeat 2 (61–93: loop 65–69, H5–H6)`, `C-terminal (94–118: H7)` — aderente al paper sperimentale di `2MYO` e rispetta l'integrità dei 2 loop $\beta$-hairpin.
  - **Su che base abbiamo scelto le 4 Candidate Regions (incrocio `RMSF` + `Q_i` + perdita `H-bonds`):**
    1. `28–37` (`LYS 37`, crollo `LEU32–ARG36` `43.6% -> 0.08%`)
    2. `46–55` (core sepolto `ILE 50, LEU 51, PHE 53` in `Repeat 1`)
    3. `60–69` (giunzione `ASP 60 – HIS 68`, crollo `THR70–ILE61` e `LYS97–ASP65`)
    4. `94–103` (ingresso `C-terminal`, crollo `THR103–ASP100` `15.6% -> 0.32%`).
- **Dove reperire dati e fonti:**
  - [`report/02_wild_type_characterization_and_baseline.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_wild_type_characterization_and_baseline.md)
  - [`SCALETTA_PROF/PUNTO_2_CARATTERIZZAZIONE_E_BASELINE_WT.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/SCALETTA_PROF/PUNTO_2_CARATTERIZZAZIONE_E_BASELINE_WT.md) & [`02_wt_baseline/RISULTATI_E_COMMENTO_FASE2.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/RISULTATI_E_COMMENTO_FASE2.md)
  - Tabelle: [`report/02_utilities/01_wt_preparation_checks.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_utilities/01_wt_preparation_checks.txt), [`report/02_utilities/07_WT_structural_region_summary.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_utilities/07_WT_structural_region_summary.csv), [`02_wt_baseline/unstable_residues_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/unstable_residues_report.txt)
- **Grafici & Comandi VMD 3D:**
  - Immagini: `report/02_utilities/02_WT_RMSD_time.png`, `03_WT_Q_time.png`, `04_WT_Rg_SASA.png`, `05_WT_Alpha_Hbonds.png`, `06_WT_Q_regions.png`, `08_WT_structural_regions_corrected.png`, `06_presentation/fig5_wt_baseline_plots.pdf`.
  - VMD 3D (box solvatato & mappa instabilità $Q_i$): `vmd -e scripts/viz_02_solvated_box.tcl` e `vmd -e scripts/viz_03_wt_instability.tcl`

---

### 3️⃣ Mutation Hypotheses and Screening (Funnel, Scarti e Risparmio Computazionale)
*Requisiti `Laboratory_Project_2026.pdf`: Slide 3 (LO2: explain how substitutions improve local interactions), Slide 4 (Methods: mutation rationale, selection), Slide 7 (computational filter su $20^{118}$), Slide 15–19 (FoldX workflow: `RepairPDB` $\rightarrow$ `PositionScan` $\rightarrow$ `BuildModel` con matched WT references $\rightarrow$ `Selection for MD`).*

- **Punti chiave da inserire:**
  - **Imbuto di riduzione dello spazio mutazionale (risparmio tempo computazionale):**
    - `118` residui totali $\rightarrow$ **-24 Functional Mask (`3AAA`)** $\rightarrow$ **-6 terminali liberi (`1–3, 116–118`)** $\rightarrow$ **`88` posizioni modificabili ($88 \times 20 = 1.760$ varianti)** $\rightarrow$ filtraggio sulle regioni instabili e consenso $\rightarrow$ `13` varianti in `FoldX BuildModel` ($n=5$) $\rightarrow$ **`4` candidati finalisti per MD**.
  - **Strumenti ortogonali utilizzati (`5 Pilastri`):**
    1. **MD WT 400 K** ($\langle Q_i \rangle$ e H-bonds locali);
    2. **FoldX 5** (`RepairPDB` clash da `69.1` a `13.5 kcal/mol` $\rightarrow$ `PositionScan` 88/88 $\rightarrow$ `BuildModel` $n=5$ con *matched WT references*);
    3. **FireProt 2.0** (ibrido Rosetta/FoldX + evolutivo *Back-to-Consensus*);
    4. **ThermoMPNN & ThermoMPNN-D** (Deep Learning GNN per singole e doppie);
    5. **Consenso Ankyrin** (*Mosavi et al., PNAS 2002* — motivo `TPLH`).
  - **Focus su COSA ABBIAMO SCARTATO e perché (cruciale per mostrare rigore scientifico!):**
    - Tutte le mutazioni nei 24 residui della **Functional Mask `3AAA`** (es. `TRP 8, ASP 44, GLU 78, HIS 67`).
    - **Falsi positivi idrofobici di superficie:** es. **`D60L`** (buon $\Delta\Delta G$ in FoldX `-2.28 kcal/mol`, ma scartato da FireProt perché sostituisce un Aspartato carico esposto sul loop con una Leucina idrofobica $\rightarrow$ rischio aggregazione).
    - **Combinazioni sub-additive in `BuildModel`:** es. `K37L + H68G` ($\Delta\Delta\Delta G_{\text{coupling}} = +1.28\text{ kcal/mol}$) scartata come doppia a favore della coppia 100% additiva `K37L + T103L` ($\Delta\Delta\Delta G = -0.06\text{ kcal/mol}$).
  - **Shortlist finale (3 mutazioni *Back-to-Consensus*, 1 per regione target):**
    - **`K37L`** (`Repeat 1`, `-2.44 ± 0.87 kcal/mol`), **`T103L`** (`C-terminal`, `-1.15 ± 0.14 kcal/mol`), **`H68G`** (`Repeat 2`, #1 BTC FireProt 78% & #1 ThermoMPNN, `-1.73 ± 0.01 kcal/mol`) $\rightarrow$ **Doppio `K37L+T103L`** (`-3.65 ± 0.66 kcal/mol`) e **Triplo `K37L+H68G+T103L`** (`-5.03 ± 0.21 kcal/mol`).
- **Dove reperire dati e fonti:**
  - [`SCALETTA_PROF/PUNTO_3_IPOTESI_MUTAZIONALI_E_SCREENING.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/SCALETTA_PROF/PUNTO_3_IPOTESI_MUTAZIONALI_E_SCREENING.md) & [`03_mutations/RISULTATI_E_COMMENTO_FASE3.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/RISULTATI_E_COMMENTO_FASE3.md)
  - [`03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md)
  - Tabelle: [`03_mutations/consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt), [`03_mutations/buildmodel_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/buildmodel_epistasis_report.txt)
- **Grafici & Comandi VMD 3D:**
  - Grafici: `06_presentation/fig6_ankyrin_consensus_alignment.pdf`, `fig7_foldx_position_scan_heatmap.pdf`, `fig8_consensus_energy_breakdown.pdf`, `fig8b_buildmodel_singles_doubles_epistasis.pdf`.
  - VMD 3D (zoom sovrapposizione WT Rosso vs Mutazioni `K37L, H68G, T103L`): `vmd -e scripts/viz_04_mutations_zoom.tcl`

---

### 4️⃣ Matched MD Comparison (Individuali $\rightarrow$ Scelta Combo Migliore $\rightarrow$ Focus `WT` vs `Triplo`)
*Requisiti `Laboratory_Project_2026.pdf`: Slide 3 (LO3: controlled computational comparison), Slide 10, 15, 21 (same protocol & multiple velocity seeds).*

- **Punti chiave da inserire:**
  - **Controllo Topologico Matched (Gate 2):** stessa box PBC (`72.59 × 58.52 × 58.88 Å`), carica netta `+0.0000 e` (`26 Na+ / 20 Cl-` per le varianti con `K37L`, `25 Na+ / 20 Cl-` per le neutre), stesso force field e stessi 3 velocity seeds (`12345, 23456, 34567`).
  - **Parte 4.1 — Risultati individuali delle 4 varianti (`K37L`, `T103L`, `Doppio`, `Triplo`) e logica di scelta:**
    - `K37L` salva il core di `Repeat 1` (`ILE 50, LEU 51, PHE 53`) ma lascia mobile `ASP 60` (`RMSF = 5.00 Å`).
    - `T103L` protegge il C-terminale e riduce RMSD e SASA.
    - `Doppio (K37L+T103L)` massimizza compattezza ($R_g$) ma conserva la flessibilità di cerniera su `ASP 60` (`RMSF = 3.35 Å`).
    - **`Triplo (K37L+H68G+T103L)` = Combo Migliore:** `H68G` rilassa il $\beta$-turn `65–69`, abbatte `RMSF(ASP 60)` a `1.50 Å` e `RMSF(LEU 37)` a `0.91 Å`, vincendo su tutte le metriche!
  - **Parte 4.2 — Focus Matched Comparison `WT` vs `Triplo` (con triplicato `n = 3` seeds $\rightarrow$ analisi di varianza $\mu \pm \sigma$):**
    - Confronto diretto delle curve $Q(t)$ (da `Q(0)=1.0`), `RMSD`, `3AAA Mask RMSD`, `RMSF` e **nuove metriche chieste dal Prof (`H-Bonds 3.5 Å / 30°` + comportamento delle 7 $\alpha$-eliche `H1–H7`, focus su `H6` e `H7`)**.
- **Dove reperire dati e fonti (aggiornati man mano da Oleg):**
  - [`04_md_simulations/matched_systems_qc_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/matched_systems_qc_report.txt)
  - [`05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md) & [`final_evaluation.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/final_evaluation.md)
  - Cartella plot individuali: [`05_analysis/single_runs/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs) (`00_WT`, `01_K37L`, `02_T103L`, `03_K37L_T103L`, `04_K37L_H68G_T103L`)
  - Cartella focus `WT vs Triplo`: [`06_presentation/presentazione_prof_WT_vs_Triplo/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/presentazione_prof_WT_vs_Triplo)
- **Grafici & Comandi VMD 3D:**
  - Grafici: `fig9_matched_q_time_and_repeats.pdf`, `fig10_per_residue_delta_qi_and_rmsf.pdf`, `01_WT_vs_Triplo_Q_value_tempo_e_moduli.png`, `02_WT_vs_Triplo_RMSD_Rgyr_SASA_3AAA.png`, `03_WT_vs_Triplo_per_residuo_Qi_e_RMSF.png`, `05_hbonds_comparison_WT_vs_Mutants.png`, `06_alpha_helices_timeline_and_H1_H7_behavior.png`.
  - VMD 3D (animazione MD 400 K del singolo sistema o sovrapposizione `WT` vs `Triplo` con eliche `H6–H7`):
    - `vmd -e scripts/viz_05_md_comparison_400K.tcl` *(opzionale: `-args WT 1` oppure `-args K37L 2`)*
    - `vmd -e scripts/viz_06_wt_vs_triple_400K.tcl` *(sovrapposizione 3D `WT` Rosso vs `Triplo` Viola a 400 K, opzionale `-args 2` o `3`)*

---

### 5️⃣ Conclusioni Finali, Incertezza Statistica & Contributi del Team
*Requisiti `Laboratory_Project_2026.pdf`: Slide 3 (LO4: assess uncertainty & discuss how computational screening supports experimental evaluation), Slide 4 (Results & interpretation + Brief description of each group member's contribution), Slide 15 (Evidence and decision).*

- **Punti chiave da inserire:**
  - **Tabella Finale Riassuntiva `WT` vs `Triplo` ($\mu \pm \sigma$ su $n=3$ seeds + $p$-value):** `Q: all-all`, `Q: hydrophilic-all`, `Q_Rep1`, `Backbone RMSD`, `Sito 3AAA RMSD`, `Rgyr`, `SASA`, `H-Bonds (3.5 Å)`, persistenza `7 α-eliche (H1–H7)` e conteggio residui critici.
  - **Interpretazione critica & limiti onesti (LO4):** valore assoluto di $Q$ vs `sdAb` di Bekker (assenza ponti disolfuro, $T_m \approx 52^\circ\text{C}$, cutoff severo `4.5 Å` su heavy atoms); scelta giustificata del protocollo single-stage NPT + finestra finale 30%.
  - **Ponte verso la validazione sperimentale (LO4):** come questo screening riduce $20^{118}$ varianti a 1–2 candidati prioritari da testare *in vitro* (es. *Circular Dichroism / DSF* per $T_m$ e saggio di binding ad *Actin-Capping Protein* per confermare la funzione).
  - **Slide finale obbligatoria (Slide 4):** breve ripartizione dei contributi di **Alessia** (caratterizzazione strutturale/SASA/H-bonds statici `2MYO`, screening Cloud FireProt 2.0 & ThermoMPNN, triplicato MD Wild-Type, redazione report e slide) e **Oleg** (setup topologico `CHARMM36m` e Functional Mask `3AAA`, screening FoldX `PositionScan` + `BuildModel`, triplicato MD mutanti, pipeline di analisi VMD/Python e confronto multi-metrica).
- **Dove reperire dati:**
  - [`05_analysis/bekker_table2_and_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/bekker_table2_and_epistasis_report.txt)
  - [`06_presentation/presentazione_prof_WT_vs_Triplo/00_TABELLA_CONFRONTO_WT_vs_TRIPLO.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/presentazione_prof_WT_vs_Triplo/00_TABELLA_CONFRONTO_WT_vs_TRIPLO.txt)
  - [`06_presentation/presentazione_prof_WT_vs_Triplo/05_hbonds_and_helices_summary.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/presentazione_prof_WT_vs_Triplo/05_hbonds_and_helices_summary.txt)
  - Grafico sintetico finale: [`06_presentation/fig11_epistatic_cycle_and_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig11_epistatic_cycle_and_observables.pdf)
