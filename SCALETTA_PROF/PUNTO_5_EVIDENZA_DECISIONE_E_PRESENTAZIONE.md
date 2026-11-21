# 📊 PUNTO 5 DELLA SCALETTA DEL PROFESSORE — "5. Evidence and decision: Compare Q and supporting observables. If necessary re-iterate"

**Riferimento Ufficiale nella Scaletta del Corso:**  
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 3** (*Learning Objective 4: Scientific interpretation — Assess uncertainty and discuss how computational screening can support experimental evaluation*), **Slide 4** (*Project outputs: Project Presentations Wednesday, December 16th, maximum 10 minutes — Introduction, Methods, Results and interpretation, Brief description of each group member's contribution*), **Slide 11–13** (*Q-value and experimental $T_m$; Higher Q supports greater contact retention under the simulated conditions*), **Slide 15** (*Project workflow — Step 5: Evidence and decision. Compare Q and supporting observables. If necessary re-iterate*), **Slide 22** (*RMSF: local flexibility -> `analysisMD.tcl`; NOTE: flexibility alone does not establish local instability or identify a beneficial mutation*), **Slide 23** (*Q-value - Native contacts*).

---

## 1. 📖 Introduzione al Problema e Riferimenti alle Slide di Teoria del Corso

Il **Punto 5** chiude il ciclo scientifico del progetto: una volta raccolte le traiettorie MD multi-seed di `WT`, `Mut1`, `Mut2` e `Double` a `400 K`, come trasformiamo i file binari `.dcd` in una **dimostrazione quantitativa, statisticamente rigorosa e biofisicamente interpretata** dell'aumento di termostabilità?

1. **L'Osservabile Primaria Definisce la Conclusione (`Laboratory_Project_2026.pdf`, Slide 11–13 & 23; *Bekker et al., 2019* pp. 431–435; *Best, Hummer & Eaton, PNAS 2013*):**
   - Nella **Slide 12 (*"The observable defines the conclusion: Higher Q supports greater contact retention under the simulated conditions"*)** e nella **Slide 13 (*"Evaluate protein variants stability by measuring the Q-value, the fraction of native contacts with respect to the wild type"*),** il Prof. Gautieri indica la **frazione di contatti atomici nativi ($Q$-value)** come metrica primaria di successo.
   - Perché il $Q$-value batte l'RMSD globale? Come dimostrato da *Bekker et al. (2019)* a pagina 431 (*"no significant difference of the average RMSDs between the antibodies is observed, suggesting that the RMSD is not a good measurement to indicate $T_m$"*), l'RMSD globale può essere dominato dallo sbandieramento di una singola coda terminale libera senza che il core della proteina si stia aprendo. Il **$Q$-value**, invece, misura direttamente la frazione di interazioni atomiche a corto raggio ($\le 4.5\text{ Å}$ tra atomi pesanti `noh`) preservate rispetto alla struttura nativa!
   - Inoltre, *Bekker et al. (2019)* (Fig. 2B–C p. 432 e Table II p. 434) dimostrano che suddividendo gli amminoacidi per classe chimico-fisica (`Lecture01`, Slide 17–20):
     - **`Q: hydrophilic-all`** (`Asp, Glu, Gln, Asn, Arg, Lys, His` vs tutti i residui) ha la **massima correlazione con la $T_m$ sperimentale ($r = 0.84$)**;
     - **`Q: all-all`** (tutti i contatti atomici pesanti) ha correlazione **$r = 0.79$**;
     - **`Q: hydrophobic-small`** (`Phe, Tyr, Trp, Leu, Val, Ile, Met, Cys, Pro` vs `Gly, Ala, Ser, Thr`) ha correlazione **$r = 0.76$**.
2. **Le Osservabili Meccaniche e Strutturali di Supporto (`Laboratory_Project_2026.pdf`, Slide 15 & 22; `Lecture04`, Slide 29–33; `analysisMD.tcl`):**
   - La **Slide 15 (*"Compare Q and supporting observables"*)** e la **Slide 22** richiedono di affiancare al $Q$-value le metriche classiche di post-processing delle traiettorie (`Lecture04`, Slide 29–33; [`laboratorio/analysisMD.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/analysisMD.tcl)):
     - **Pre-processing obbligatorio (`Slide 22`):** ricomposizione delle molecole periodiche (`pbc wrap`) e allineamento roto-traslazionale sul backbone del frame iniziale (`measure fit`);
     - **RMSD del backbone (`Lecture04`, Slide 30):** spostamento quadratico medio nel tempo;
     - **RMSF dei $C_\alpha$ per residuo (`Lecture04`, Slide 31; `Slide 22`):** fluttuazione locale attorno alla posizione media (ricordando il monito della Slide 22: *"flexibility alone does not establish local instability or identify a beneficial mutation"*, motivo per cui confrontiamo $\Delta\text{RMSF} = \text{RMSF}_{\text{Mut}} - \text{RMSF}_{\text{WT}}$ insieme a $\Delta\langle Q_i \rangle$!);
     - **Struttura Secondaria nel tempo (`Lecture04`, Slide 32; `analysisMD.tcl`):** frazione di residui in $\alpha$-elica calcolata tramite STRIDE (`mol ssrecalc`);
     - **Legami a Idrogeno e Ponti Salini (`Lecture04`, Slide 33; `analysisMD.tcl`):** numero totale di H-bond nel tempo e occupancy (%) delle singole coppie donatore–accettore;
     - **Raggio di Girazione ($R_g$) e SASA (`analysisMD.tcl`):** misura della compattezza globale e dell'esposizione al solvente.
3. **Valutazione dell'Incertezza Statistica e Rilevanza Sperimentale (`Laboratory_Project_2026.pdf`, Slide 3 — *Learning Objective 4*):**
   - Ogni metrica va espressa come **Media $\pm$ Deviazione Standard ($\mu \pm \sigma$)** calcolata tra le $n = 3$ repliche indipendenti (*velocity seeds*) nella finestra finale equilibrata, discutendo se l'effetto del Doppio Mutante sia **additivo** o **sinergico/epistatico** (*Bekker et al., 2019* p. 434) e come lo screening *in silico* guidi la successiva validazione sperimentale *in vitro* (es. DSF / Dicroismo Circolare sulla curva di thermal unfolding della Slide 5!).

---

## 2. 🎯 Cosa Cerchiamo nel Punto 5 (Gli Obiettivi dell'Analisi Finale)

1. **Dimostrazione Primaria (*Bekker et al. Table II*):** Dimostrare che `Mut1`, `Mut2` e soprattutto `Double` mantengono a `400 K` valori medi di **`Q: hydrophilic-all`**, **`Q: all-all`** e **`Q: hydrophobic-small`** superiori al Wild-Type, con riduzione della deviazione standard tra le repliche e riduzione del numero totale di residui localmente instabili ($\langle Q_i \rangle_{400\text{K}} < 0.60$).
2. **Dimostrazione Biofisica Specifica per la Myotrophin (*Lowe & Itzhaki, PNAS 2007*):**
   - Misurare il guadagno di contatti nativi nei singoli Ankyrin Repeats (**`Q_ANK1`**, **`Q_ANK2`**, **`Q_ANK3`**, **`Q_ANK4`**) e nelle **3 interfacce inter-repeat (`Q_ANK1-2`, `Q_ANK2-3`, `Q_ANK3-4`)**;
   - Verificare se la combinazione delle mutazioni nel Doppio Mutante (es. `K37L` in `ANK2` + `T103L` in `ANK4`) stabilizza per **propagazione cooperativa allosterica** anche il repeat centrale `ANK3`!
3. **Conferma Meccanica e Chimica (`analysisMD.tcl`):**
   - Mostrare l'abbattimento dei picchi di $\text{RMSF}_{C_\alpha}$ nei loop/eliche mutati ($\Delta\text{RMSF} < 0$), la maggiore tenuta della `% Alpha-Helix` e del raggio di girazione $R_g$, l'eliminazione del ponte salino parassita `LYS 37 -> ASP 60` e la preservazione dell'integrità geometrica dei **24 residui della Functional Mask (`3AAA`)**.
4. **Sintesi Esecutiva per la Presentazione del 16 Dicembre (`Slide 4`, max 10 minuti):** Impaginare risultati, figure vettoriali `.pdf`, render 3D VMD e divisione chiara dei contributi di Oleg e Alessia.

---

## 3. 🛠️ Come Lo Cercheremo (La Politica di Validazione a 5 Livelli Già Programmata)

L'intero protocollo matematico e statistico del Punto 5 è stato progettato passo-passo in **[`ROADMAP/ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md)** e produrrà nella cartella `05_analysis/` e `06_presentation/`:

1. **Tabella Master Comparativa (`05_analysis/bekker_table2_comparison.txt` & `.csv`):**
   - Replica esatta della **Table II di Bekker et al. (2019, p. 434)** estesa con `Q_ANK1..4`, `% Alpha-Helix`, $R_g$, numero di residui instabili ($Q_i < 0.60$) e $p$-value del test $t$ di Welch rispetto al WT.
2. **4 Figure Vettoriali Finali in `.pdf` per le Slide di Risultati (`06_presentation/`):**
   - **`fig9_matched_q_time_and_repeats.pdf`:** Curve temporali $Q(t)$ a `400 K` con bande di confidenza $\mu \pm \sigma$ sui 3 seed per `WT` (rosso), `Mut1` (blu), `Mut2` (arancio) e `Double` (verde), più barplot con barre d'errore di `Q: hydrophilic-all`, `Q: all-all`, `Q: hydrophobic-small` e `Q_ANK1..4`.
   - **`fig10_per_residue_delta_qi_and_rmsf.pdf`:** Profili sovrapposti lungo i 118 residui di $\langle Q_i \rangle_{400\text{K}}$ e $\Delta\text{RMSF}_{C_\alpha}$, mostrando visivamente la "guarigione" dei residui instabili sotto la soglia `0.60`.
   - **`fig11_epistatic_cycle_and_contact_map.pdf`:** Ciclo di additività/epistasi del Doppio Mutante ($\Delta Q_{1+2}$ vs $\Delta Q_1 + \Delta Q_2$) e heatmap 2D della matrice di differenza dei contatti $\Delta C_{ij} = \langle C_{ij}^{\text{Double}} \rangle - \langle C_{ij}^{\text{WT}} \rangle$.
   - **`fig12_analysisMD_supporting_observables.pdf`:** Confronto multi-seed delle osservabili di `analysisMD.tcl` (`% Alpha-Helix`, `Hydrogen Bonds`, `Rgyr`, `SASA` e `RMSD della Functional Mask 3AAA`).

---

## 4. 💡 Cosa Abbiamo Già Trovato e Predisposto Finora

- **Baseline WT a `300 K` e `400 K` già completata e analizzata (`02_wt_baseline/`):**
  - Tutti i valori di riferimento del Wild-Type per la finestra equilibrata `876–1251` sono già calcolati: `RMSD = 7.10 Å`, `Q_all = 0.693`, `Q_phil = 0.678`, `Q_phob_small = 0.625`, `Q_ANK1 = 0.703`, `Q_ANK2 = 0.594`, `Q_ANK3 = 0.662`, `Q_ANK4 = 0.675`, e **`28` residui dello scaffold interno con $\langle Q_i \rangle_{400\text{K}} < 0.60$**.
- **5 Figure Vettoriali `.pdf` su 9 già pronte in [`06_presentation/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation):**
  1. [`fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf) (Caratterizzazione WT: RMSD, $Q(t)$, RMSF, $\langle Q_i \rangle_{400\text{K}}$).
  2. [`fig5b_wt_analysisMD_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5b_wt_analysisMD_observables.pdf) (Osservabili `analysisMD.tcl` del WT: `% Alpha-Helix`, H-bonds, $R_g$, $Q$ per Repeat).
  3. [`fig6_ankyrin_consensus_alignment.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig6_ankyrin_consensus_alignment.pdf) (Allineamento dei 4 Ankyrin Repeats al consenso di *Mosavi et al., 2002* con Functional Mask `3AAA`).
  4. [`fig7_foldx_position_scan_heatmap.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig7_foldx_position_scan_heatmap.pdf) (Top 20 mutazioni stabilizzanti da FoldX `PositionScan` `88/88`).
  5. [`fig8_consensus_energy_breakdown.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8_consensus_energy_breakdown.pdf) (Scomposizione fisico-chimica dei migliori candidati a 5 pilastri).

---

## 5. 🎤 Architettura della Presentazione Finale del 16 Dicembre (`Slide 4`: Max 10 Minuti, 10 Slide)

La **Slide 4 di [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf)** richiede esplicitamente 4 elementi in massimo **10 minuti**:
1. *Introduction on the target protein*
2. *Methods (mutation rationale, selection, ...)*
3. *Results and interpretation*
4. *Brief description of each group member's contribution*

Ecco la struttura a **10 Slide (~1 minuto a slide, 5 minuti Oleg + 5 minuti Alessia)** già allineata ai 4 criteri di valutazione del Prof. Gautieri (*scientific reasoning, methodological care and communication*):

| Slide | Titolo & Contenuto Scientifico | Figure / Tabelle Pronte o Previste | Chi Parla |
| :---: | :--- | :--- | :---: |
| **1** | **Titolo, Team #10 & Divisione dei Contributi (`Slide 4`)**<br>Panoramica del workflow a 5 punti e ripartizione paritaria: **Oleg** (Setup CHARMM36m, Benchmark HPC `+p8` AVX-512, MD WT/Mutanti, FoldX 88/88, Analisi VMD/$Q$-value) & **Alessia** (Analisi strutturale `3AAA` & Consenso Ankyrin, Screening Cloud FireProt 2.0, ThermoMPNN-D Singole/Doppie, Validazione Epistatica e Interpretazione Biofisica). | Schema grafico del workflow + Tabella ruoli | **Oleg & Alessia** (30s) |
| **2** | **Punto 1A — Identikit Biofisico di Myotrophin (`2MYO`) e Complesso Funzionale (`3AAA`)**<br>Architettura a 4 Ankyrin Repeats (`ANK1–4`), confronto Ratto NMR (`2MYO`) vs Uomo X-Ray (`3AAA`, $99.2\%$ identità, $\text{RMSD} = 2.365\text{ Å}$) e mappatura oggettiva della **Functional Mask (24 residui a $< 4.5\text{ Å}$ da Actin-Capping Protein)**. | Render 3D VMD da [`viz_01_structure_and_mask.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_01_structure_and_mask.tcl) | **Alessia** (1m) |
| **3** | **Punto 1B & 2A — Cura Metodologica: Rete H-Bond delle Istidine, Setup `CHARMM36m` & Benchmark HPC**<br>Scoperta chimica su `HIS 40 (HSD -> SER 74)` e `HIS 67 (HSE -> GLY 35 / LYS 37)`, box `TIP3P` (`21.461` atomi, $0.15\text{ M NaCl}$), speedup **`+40%` su 8 core fisici (`12.68 ns/day`)** e tabella del **Computational Budget (`Slide 21`)**. | Render 3D box [`viz_02_solvated_box.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_02_solvated_box.tcl) + Tabella Budget | **Oleg** (1m) |
| **4** | **Punto 2B — Caratterizzazione Dinamica del Wild-Type (`300 K` Controllo vs `400 K` Stress Termico)**<br>Validazione a `300 K` ($R_g = 14.5\text{ Å}, \alpha\text{-elica } 46\%$) e diagnosi del cedimento a `400 K`: `ANK2` è il modulo più fragile ($Q_{\text{ANK2}} = 0.586$), collasso del ponte `LYS 37 / ASP 60` ($\langle Q_{37} \rangle = 0.485$) e dell'ingresso C-cap `ANK4` (`LEU 102 – THR 103`, $\langle Q_{103} \rangle = 0.508$). | [`fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf) + [`fig5b_wt_analysisMD_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5b_wt_analysisMD_observables.pdf) + Render 3D [`viz_03_wt_instability.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_03_wt_instability.tcl) | **Oleg** (1.5m) |
| **5** | **Punto 3A — "Orthogonal Consensus Funnel": Consenso Ankyrin, FireProt 2.0 & ThermoMPNN-D**<br>Allineamento ai 33 residui canonici (*Mosavi et al., 2002*), screening evolutivo *Back-to-Consensus* di FireProt 2.0 (`H68G` #1 BTC) e Graph Neural Network ThermoMPNN-D (`K37` #1 assoluto). Perché abbiamo scartato la "trappola" `S74A` (protezione ponte H `HIS 40 -> SER 74`) e gli *helix-breakers* `E52P`/`T103G`. | [`fig6_ankyrin_consensus_alignment.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig6_ankyrin_consensus_alignment.pdf) + Tabella Filtri | **Alessia** (1.5m) |
| **6** | **Punto 3B — Screening FoldX (`88/88`), Energy Breakdown & Modelli 3D (`BuildModel`)**<br>Risultati di `RepairPDB` ($\Delta E = -149.4\text{ kcal/mol}$), `PositionScan` sul 100% dello scaffold e scomposizione energetica dei candidati a **5/5 pilastri** (**`K37L`** in `ANK2`, **`T103L`** in `ANK4`, **`H68G`** in `ANK3`) + zoom 3D di `BuildModel`. | [`fig7_foldx_position_scan_heatmap.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig7_foldx_position_scan_heatmap.pdf) + [`fig8_consensus_energy_breakdown.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8_consensus_energy_breakdown.pdf) + Zoom 3D VMD | **Alessia & Oleg** (1m) |
| **7** | **Punto 4 & 5A — Matched MD Comparison (`400 K`, $n=3$ Velocity Seeds): Risultati $Q$-Value (*Bekker Table II*)**<br>Curve $Q(t)$ multi-seed ($\mu \pm \sigma$) e **Table II** completa (`Q: hydrophilic-all`, `Q: all-all`, `Q: hydrophobic-small` e riduzione del numero di residui instabili $\langle Q_i \rangle < 0.60$). | `fig9_matched_q_time_and_repeats.pdf` + **Table II** | **Oleg** (1m) |
| **8** | **Punto 5B — Stabilità Risolta per Ankyrin Repeat (`Q_ANK1..4`), Ciclo Epistatico & Mappe $\Delta C_{ij}$**<br>Propagazione cooperativa della stabilità dal nucleo di folding C-cap (`ANK4`, *Lowe & Itzhaki 2007*) lungo il solenoide e test di additività/sinergia del Doppio Mutante (`K37L + T103L`). | `fig10_per_residue_delta_qi_and_rmsf.pdf` + `fig11_epistatic_cycle_and_contact_map.pdf` | **Alessia** (1m) |
| **9** | **Punto 5C — Osservabili Meccaniche (`analysisMD.tcl`) & Preservazione del Sito Funzionale (`3AAA`)**<br>Riduzione di $\Delta\text{RMSF}_{C_\alpha}$, conservazione della `% Alpha-Helix` e di $R_g$, occupancy H-bond e verifica che l'interfaccia di legame per l'Actin-Capping Protein rimanga intatta. | `fig12_analysisMD_supporting_observables.pdf` | **Oleg** (1m) |
| **10** | **Conclusioni Scientifiche, Limiti del Modello & Prospettive Sperimentali (*In Vitro* DSF/CD)**<br>Sintesi del guadagno di stabilità ottenuto, discussione critica dell'incertezza statistica (*Learning Objective 4*) e proposta di validazione sperimentale tramite *Differential Scanning Fluorimetry (DSF)* / *Circular Dichroism*. | Box riassuntivo & Bibliografia essenziale | **Alessia & Oleg** (30s) |
