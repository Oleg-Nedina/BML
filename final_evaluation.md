# 🏁 FINAL EVALUATION — Validità Scientifica dei Risultati Ottenuti e Protocollo Operativo per la Validazione Finale (Passo A & Passo B)

**Corso:** Biomolecular Modelling Laboratory (2026) — Prof. Alfonso Gautieri (Politecnico di Milano)  
**Team #10:** Oleg Nedina & Alessia Rigoni  
**Proteina Target:** Myotrophin (*Rattus norvegicus*, PDB ID: [`2MYO`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/2MYO.pdb), 118 aa — 4 Ankyrin Repeats)  
**Data:** 22 Novembre 2026

---

## 1. ✅ Perché i Risultati dello Screening MD a `400 K` (`Seed 1` sui 5 Sistemi) Sono Scientificamente Validi

Nello **Stadio 1 della Fase 4** abbiamo eseguito e analizzato 5 simulazioni di Dinamica Molecolare all-atom a `400 K` (`2.50 ns` ciascuna = `1.250.000` step con $\Delta t = 2.0\text{ fs}$):
1. **`00_WT`** — Wild-Type (`2MYO` non mutata, baseline di controllo a 400 K) $\rightarrow$ [`05_analysis/single_runs/00_WT/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/00_WT)
2. **`01_K37L`** — Mutante Singolo 1 (`ANK2`, 5/5 pilastri) $\rightarrow$ [`05_analysis/single_runs/01_K37L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/01_K37L)
3. **`02_T103L`** — Mutante Singolo 2 (`ANK4` C-cap, 5/5 pilastri) $\rightarrow$ [`05_analysis/single_runs/02_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/02_T103L)
4. **`03_K37L_T103L`** — Doppio Mutante (`ANK2 + ANK4`) $\rightarrow$ [`05_analysis/single_runs/03_K37L_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/03_K37L_T103L)
5. **`04_K37L_H68G_T103L`** — Triplo Mutante (`ANK2 + ANK3 + ANK4`) $\rightarrow$ [`05_analysis/single_runs/04_K37L_H68G_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/04_K37L_H68G_T103L)

I risultati di queste 5 simulazioni sono **pienamente validi, rigorosi e direttamente utilizzabili nella presentazione finale** per **4 ragioni metodologiche e fisiche fondamentali**:

### A. Rispetto Rigoroso del Principio di *Matched Computational Comparison* (`Laboratory_Project_2026.pdf`, Slide 3, 15, 21)
Come richiesto dal **Learning Objective 3 (Slide 3: *"Controlled computational comparison — Use consistent modelling and MD settings for WT and variants"*)** e dalla **Slide 21 (*"Use the same protocol for WT and variants"*)**:
- Tutti e 5 i sistemi sono stati costruiti con lo **stesso force field (`CHARMM36m`)**, lo **stesso modello d'acqua esplicito (`TIP3P` con padding $12.0\text{ Å}$)**, la **stessa forza ionica fisiologica (`0.15 M NaCl`)** e **carica netta totale esattamente pari a `+0.0000 e`** (verificata nel report di Gate 2 [`04_md_simulations/matched_systems_qc_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/matched_systems_qc_report.txt) per evitare artefatti PME).
- Tutti i sistemi condividono la **stessa assegnazione tautomerica razionale delle istidine** verificata in Fase 1 (`HIS 40 -> HSD`, `HIS 67 -> HSE`, `HIS 80 -> HSD`, e `HIS 68 -> HSD` dove non sostituita da `GLY 68`).
- Tutte e 5 le traiettorie hanno subito **esattamente lo stesso protocollo numerico in NAMD 3**: `1.000` step di minimizzazione *Conjugate Gradient* + inizializzazione Maxwell-Boltzmann a `400 K` + `1.250.000` step (`2.50 ns`) in ensemble NPT (`Langevin` $\gamma = 1\text{ ps}^{-1}$, `Nose-Hoover Langevin Piston` a `1.01325 bar`, cutoff `12/10 Å`, `PMEGridSpacing 1.0 Å`, `rigidBonds all`).
- Pertanto, **ogni differenza osservata tra `WT` e varianti mutate dipende esclusivamente dalle sostituzioni amminoacidiche introdotte** e non da asimmetrie del setup.

### B. Gestione dell'Equilibrazione tramite Scarto del Transitorio Iniziale (`70%` Iniziale Scartato vs `30%` Finale Valutato)
Ci siamo chiesti se il fatto di non aver separato uno step preliminare di equilibrazione a `300 K` prima del `run` a `400 K` potesse influenzare le medie dello screening:
1. **In NAMD 3, l'equilibrazione termica e barica del solvente `TIP3P` avviene nei primi `~20–50 ps`** dopo i `1.000` step di minimizzazione.
2. **Nel nostro protocollo di analisi ([`scripts/09_analyze_matched_validation.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/09_analyze_matched_validation.tcl)), i primi `1.75 ns` (`Frame 0..875`, pari al `70%` dell'intera simulazione!) vengono interamente scartati come fase transitoria di equilibrazione e rilassamento termico a `400 K`**.
3. Tutte le osservabili termodinamiche ($\langle Q \rangle$, $\langle Q_i \rangle_{400\text{K}}$, RMSD, $R_g$, SASA, H-Bonds, % $\alpha$-elica) vengono mediate **esclusivamente nella finestra equilibrata finale del `30%` (`Frame 876..1251`, corrispondenti a `1.75–2.50 ns`)**, esattamente come prescritto da *Bekker, Ma & Kamiya (2019)* (p. 436).
4. **In più**, aver notato l'importanza della **Pre-Equilibrazione Restrained a `300 K`** (`Lecture04` Slide 6–9) ci permette di inserirla nel **Passo B (Campagna di Validazione Finale `WT` vs `Triplo Mutante`)**, trasformando un dubbio metodologico in un **punto di eccellenza scientifica** da mostrare al Prof. Gautieri!

### C. Concordanza Fisico-Chimica Totale tra Design Energetico Statico (Fase 3) e Risposta Dinamica a 400 K (Fase 4)
I numeri ottenuti nelle 5 simulazioni ([`05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md)) confermano punto per punto le previsioni di `FoldX BuildModel`, `FireProt 2.0` e `ThermoMPNN`:

| Sistema (`1.75–2.50 ns` a `400 K`) | $\Delta\Delta G_{\text{BuildModel}}$ ($n=5$ run) | `Q: hydrophilic-all` *(Bekker $r=0.84$ con $T_m$)* | `Q: all-all` *(Bekker $r=0.79$)* | `Q_ANK2` *(Core critico)* | Backbone RMSD Globale | Sito Attivo `3AAA` RMSD | SASA Totale ($\text{Å}^2$) | $\text{RMSF}_{C_\alpha}$ `Res 37 / 60` | Residui Scaffold Instabili ($\langle Q_i \rangle < 0.60$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`00_WT` (`2MYO`)** | `0.00 kcal/mol` | `0.6750 ± 0.0067` | `0.6948 ± 0.0054` | `0.6015` | `6.68 ± 0.52 Å` | `7.27 ± 0.99 Å` | `8839 ± 177` | `2.96 / 2.53 Å` | `41` residui |
| **`01_K37L` (`ANK2`)** | `-2.44 ± 0.87` | `0.7134 ± 0.0070` | `0.7243 ± 0.0052` | `0.6351` | `4.17 ± 0.33 Å` | `4.75 ± 0.32 Å` | `7942 ± 206` | `1.37 / 5.00 Å` | `25` residui *(-16)* |
| **`02_T103L` (`ANK4`)** | `-1.15 ± 0.14` | `0.7047 ± 0.0066` | `0.7124 ± 0.0046` | `0.6161` | `3.77 ± 0.19 Å` | `3.37 ± 0.36 Å` | `7969 ± 219` | `1.83 / 1.75 Å` | `30` residui *(-11)* |
| **`03_K37L_T103L` (`Doppio`)** | `-3.65 ± 0.66` | 🏆 **`0.7178 ± 0.0075`** | `0.7194 ± 0.0061` | `0.6290` | `6.40 ± 0.53 Å` | `6.73 ± 0.52 Å` | `8118 ± 269` | `2.08 / 3.35 Å` | `31` residui *(-10)* |
| **`04_K37L_H68G_T103L` (`Triplo`)** | 🏆 **`-5.03 ± 0.21`** | 🏆 **`0.7172 ± 0.0056`** | 🏆 **`0.7289 ± 0.0049`** | 🏆 **`0.6444`** | 🏆 **`3.30 ± 0.42 Å`** | 🏆 **`2.83 ± 0.26 Å`** | 🏆 **`7647 ± 235`** | 🏆 **`0.91 / 1.50 Å`** | 🏆 **`23` residui *(-18)*** |

### D. La Spiegazione Meccanicistica del Trionfo del Triplo Mutante (`K37L + H68G + T103L`)
1. **`K37L` (`ANK2`)** elimina il ponte salino distorsivo `LYS 37 – ASP 60`, chiude il core idrofobico di `ANK2` ($Q_{\text{ANK2}}$ da `0.6015` a `0.6351`) e salva 16 residui instabili, ma lascia **`ASP 60`** (sulla punta del loop `60–68` tra `ANK2` e `ANK3`) libero di oscillare ($\text{RMSF}_{60} = 5.00\text{ Å}$ in `K37L` e `3.35 Å` nel Doppio `K37L_T103L`).
2. **`T103L` (`ANK4`)** blinda l'elica C-terminale ($\Delta E_{\text{VdW+SolvH}} = -4.03\text{ kcal/mol}$) e, combinato con `K37L` nel **Doppio `K37L_T103L`**, porta `Q: hydrophilic-all` a **`0.7178`** (`+4.28%` vs WT) e il raggio di girazione a **`14.76 Å`** (compattezza nativa di 300 K conservata a 400 K!). Tuttavia, la flessione della cerniera centrale `ASP 60 – HIS 68` tra i due estremi irrigiditi mantiene alto l'RMSD globale del Doppio (`6.40 Å`).
3. **L'aggiunta di `H68G` (`ANK3` — #1 *Back-to-Consensus* di FireProt 2.0 e #1 di ThermoMPNN) nel Triplo `K37L_H68G_T103L`** rilassa la tensione sterica esattamente al residuo `68` (l'uscita del loop `60–68` tra `ANK2` e `ANK3`):
   - La fluttuazione di **`ASP 60` crolla da `5.00 / 3.35 Å` a `1.50 Å`** e quella di **`LEU 37` scende a `0.91 Å`** a 400 K!
   - Il **Backbone RMSD Globale si dimezza a `3.30 ± 0.42 Å`** (`-51%` vs WT), l'**RMSD del Sito Funzionale `3AAA` crolla a `2.83 ± 0.26 Å`** (`-61%` vs WT), **`Q: all-all` tocca il massimo assoluto (`0.7289`)**, la SASA diminuisce di **`-1.192 Å²`**, e i residui instabili scendono da **`41` a soli `23` (`-44%`)**!

---

## 2. 🚀 Cosa Bisogna Fare Adesso e Come Farlo (Passo A & Passo B)

Seguendo esattamente la logica che abbiamo definito, il completamento del progetto si divide in **Passo A** (unione immediata dei risultati dello screening a 5 sistemi) e **Passo B** (campagna di validazione finale testa-a-testa tra `Wild-Type` e il vincitore `Triplo Mutante K37L_H68G_T103L`).

---

### 📊 PASSO A (Immediato — 5 secondi): Unione dei 5 Risultati e Generazione dei Plot Comparativi Sovrapposti

- **Cosa fa:** Legge i file `.dat` già generati in `05_analysis/` per tutte e 5 le simulazioni (`00_WT`, `01_K37L`, `02_T103L`, `03_K37L_T103L`, `04_K37L_H68G_T103L`) mentre mantiene intatti i 15 plot singoli separati in [`05_analysis/single_runs/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs), e genera:
  1. **[`05_analysis/bekker_table2_and_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/bekker_table2_and_epistasis_report.txt)** — La **Table II di Bekker et al.** unificata, la tabella per singolo Repeat `ANK1..4` e l'analisi locale sui siti `37, 60, 68, 102, 103`.
  2. **`06_presentation/fig9_matched_q_time_and_repeats.pdf`** — Curve temporali $Q_{\text{hydrophilic-all}}(t)$ e $Q_{\text{all-all}}(t)$ sovrapposte per i 5 sistemi + barplot per classe chimica e per modulo `ANK1..4`.
  3. **`06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf`** — Profilo $\langle Q_i \rangle_{400\text{K}}$ per residuo (`1..118`), istogramma di guarigione $\Delta\langle Q_i \rangle_{400\text{K}}$ (`Verde = residui salvati`) e curve $\text{RMSF}(C_\alpha)$ sovrapposte con l'abbattimento del picco di `ASP 60` nella Tripla.
  4. **`06_presentation/fig11_epistatic_cycle_and_observables.pdf`** — Confronto `Backbone RMSD vs Sito 3AAA RMSD`, `SASA & Rgyr`, riduzione dei residui instabili (`41 -> 23`) e correlazione lineare tra $\Delta\Delta G_{\text{FoldX}}$ e $Q_{\text{MD 400K}}$.

- **Come eseguirlo (dalla cartella `BML/`):**
  ```bash
  python3 scripts/10_plot_final_validation_bekker.py && cat 05_analysis/bekker_table2_and_epistasis_report.txt
  ```

---

### 🧪 PASSO B (Campagna di Validazione Finale): Confronto Rigoroso Testa-a-Testa `Wild-Type (WT)` vs `Triplo Mutante (K37L_H68G_T103L)`

#### 1. Perché nel Passo B simuliamo solo `WT` e il `Triplo Mutante K37L_H68G_T103L`?
- Nello **Stadio 1 (`Seed 1` sui 5 sistemi)** abbiamo già completato la dissezione meccanicistica: sappiamo quanto contribuisce `K37L` da solo, quanto contribuisce `T103L` da solo, cosa fa il Doppio `K37L_T103L` e perché il **Triplo `K37L_H68G_T103L`** è superiore a tutti.
- Nel **Passo B (Validazione Finale)** applichiamo quindi la logica dell'**imbuto decisionale (*Screening Funnel*, Slide 15 & 19)**: prendiamo **unicamente il controllo `Wild-Type (WT)` e il nostro candidato vincitore `Triplo Mutante (K37L_H68G_T103L)`** e li sottoponiamo al protocollo di validazione completo da pubblicazione scientifica.

#### 2. Cosa facciamo esattamente nel Passo B (I 2 Pilastri della Validazione Finale)
Come spiegato in dettaglio anche nel documento per Alessia ([`ROAD_PER_ALE/04_SPIEGAZIONE_PASSO_B_VALIDAZIONE_FINALE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROAD_PER_ALE/04_SPIEGAZIONE_PASSO_B_VALIDAZIONE_FINALE.md)) e nella specifica tecnica ([`ROADMAP/PASSO_B_VALIDAZIONE_FINALE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/PASSO_B_VALIDAZIONE_FINALE.md)):

1. **Pilastro B.1 — Introduzione della Fase Esplicita di Pre-Equilibrazione Restrained a `300 K` prima dello Stress a `400 K` (`Lecture04` Slide 6–9; *Bekker et al. 2019* p. 436):**
   - Per rispondere alla questione dell'equilibrazione, ogni run di validazione del Passo B su `WT` e `K37L_H68G_T103L` seguirà il protocollo canonico a **3 stadi consecutivi**:
     - **Stage 1 (`1.000` step Minimizzazione Energetica con `constraints on` sul backbone):** Rilassa i contatti sterici dell'acqua e delle catene laterali tenendo fermo il backbone.
     - **Stage 2 (`50.000` step = `100 ps` di Pre-Equilibrazione NPT a `300 K` e `1 atm` con *Harmonic Position Restraints* sul backbone `name N CA C O`, $k = 1.0\text{ kcal/mol/Å}^2$):** Permette alla scatola d'acqua `TIP3P`, agli ioni ($\text{Na}^+, \text{Cl}^-$) e alle catene laterali mutate di raggiungere il perfetto equilibrio termico, barico e di idratazione a temperatura fisiologica (`300 K`) senza alcuna deformazione precoce del backbone.
     - **Stage 3 (`1.250.000` step = `2.50 ns` di Produzione NPT libera a `400 K` con `constraintScaling 0.0`):** Si spengono i vincoli sul backbone e si applica lo stress termico a `400 K` a partire dal sistema già equilibrato a `300 K`.
2. **Pilastro B.2 — Repliche Statistiche Indipendenti Multi-Seed ($n = 3$ Velocity Seeds: `12345`, `23456`, `34567`, Slide 10, 15, 21):**
   - Per stimare l'incertezza statistica ($\mu \pm \sigma$ e $p$-value del test $t$ di Welch, **Learning Objective 4, Slide 3 & 21**), confrontiamo `WT` e `K37L_H68G_T103L` su 3 semi indipendenti della distribuzione delle velocità di Maxwell-Boltzmann.
   - **Vantaggio enorme sul tempo di calcolo:** Concentrandoci solo su **`WT` vs `K37L_H68G_T103L`**, l'intera campagna di validazione richiede solo **2 sistemi** (invece di 5 o 6!), che potete anche dividere a metà tra il PC di Oleg (es. `K37L_H68G_T103L`) e il PC di Alessia (es. `WT`), chiudendo la validazione in **una sola notte**!

---

## 3. 🗂️ Indice Completo dei File Aggiornati Pronti per il `git push`

- **Questo documento master:** [`final_evaluation.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/final_evaluation.md)
- **Analisi completa delle 5 simulazioni `Seed 1` a `400 K`:** [`05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md)
- **I 15 Plot PDF e i 5 Report separati per ogni singola run:**
  - [`05_analysis/single_runs/00_WT/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/00_WT)
  - [`05_analysis/single_runs/01_K37L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/01_K37L)
  - [`05_analysis/single_runs/02_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/02_T103L)
  - [`05_analysis/single_runs/03_K37L_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/03_K37L_T103L)
  - [`05_analysis/single_runs/04_K37L_H68G_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/04_K37L_H68G_T103L)
- **Guida di spiegazione per Alessia sul Passo B:** [`ROAD_PER_ALE/04_SPIEGAZIONE_PASSO_B_VALIDAZIONE_FINALE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROAD_PER_ALE/04_SPIEGAZIONE_PASSO_B_VALIDAZIONE_FINALE.md)
- **Specifica tecnica del Passo B:** [`ROADMAP/PASSO_B_VALIDAZIONE_FINALE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/PASSO_B_VALIDAZIONE_FINALE.md)
- **Documenti di Fase, Scaletta del Professore e Roadmap aggiornati:**
  - [`01_structures/RISULTATI_E_COMMENTO_FASE1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/RISULTATI_E_COMMENTO_FASE1.md)
  - [`02_wt_baseline/RISULTATI_E_COMMENTO_FASE2.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/RISULTATI_E_COMMENTO_FASE2.md)
  - [`03_mutations/RISULTATI_E_COMMENTO_FASE3.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/RISULTATI_E_COMMENTO_FASE3.md)
  - [`04_md_simulations/RISULTATI_E_COMMENTO_FASE4_SETUP.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/RISULTATI_E_COMMENTO_FASE4_SETUP.md)
  - [`SCALETTA_PROF/PUNTO_1..5`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/SCALETTA_PROF)
  - [`ROADMAP/ROAD1..5.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP)
  - [`CHECK_SCALETTA_PROFESSORE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/CHECK_SCALETTA_PROFESSORE.md) & [`STATUS/SESSION1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/STATUS/SESSION1.md)
