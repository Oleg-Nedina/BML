# 📋 SCALETTA.md — Punto della Situazione, Logica Scientifica e Piano Operativo Completo

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Proteina Target:** Myotrophin (*Rattus norvegicus*, PDB ID: `2MYO`, 118 amminoacidi)  
**Presentazione Finale:** Mercoledì 16 Dicembre (max 10 minuti)  
**Riferimenti Ufficiali:** [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf), *Bekker, Ma & Kamiya (2019)* (`Protein Science 28:429–438`), *Mosavi et al. (2002)* (`PNAS 99:16029–16034`), *Takeda et al. (2010)* (`PLoS Biology 8:e1000416`).

---

## 1. 🔍 Punto della Situazione Attuale (Dove Siamo Davvero e Perché C'era Confusione)

### Perché ci siamo sentiti persi?
Il flusso logico del progetto prevede **5 Fasi sequenziali** (`Fase 1` $\rightarrow$ `Fase 5`).  
Siccome la simulazione vera della **Fase 2** (Dinamica Molecolare del Wild-Type a `300 K` e `400 K`) dura circa **7.5 ore** ed è stata lasciata per la notte, abbiamo **anticipato i calcoli statici della Fase 3** (`FoldX RepairPDB`, `PositionScan`, e `Consenso Ankyrin`) che non dipendono dalla CPU di NAMD.  
Questo "salto temporale" — insieme agli oltre **1.000 file `.pdb` temporanei** sputati da FoldX dentro `03_mutations/` — ha reso difficile vedere a colpo d'occhio cosa è già concluso, cosa era solo un test e cosa manca per andare avanti.

Ecco lo **stato reale verificato** di ogni cartella del progetto:

| Cartella | Fase Corrispondente | Stato Reale Attuale | Dettaglio Preciso (Zero Ipotesi) |
| :--- | :--- | :--- | :--- |
| **`01_structures/`** | **Fase 1** (`ROAD1.md`): Studio Proteina & Setup Scatola d'Acqua | ✅ **COMPLETATA AL 100%** | • Scaricati `2MYO.pdb` (Ratto NMR, 118 aa) e `3AAA_complex_ref.pdb` (complesso funzionale).<br>• Verificate 3 Cisteine (0 ponti disolfuro) e assegnate le 4 Istidine (`HIS 40 -> HSD`, `HIS 67 -> HSE`, `HIS 68/80 -> HSD`).<br>• Estratta la **Functional Mask** oggettiva: **24 residui intoccabili** ($< 4.5\text{ Å}$ da Actin-Capping Protein in `3AAA`).<br>• Costruito sistema in **CHARMM36m + acqua `TIP3P` ($12\text{ Å}$) + $0.15\text{ M NaCl}$** (`2myo_wb_ionized.psf/.pdb`, **21.461 atomi**).<br>• Foto VMD della Fase 1 già scattate. |
| **`02_wt_baseline/`** | **Fase 2** (`ROAD2.md`): Benchmark & MD Baseline del Wild-Type | 🔄 **SETUP E TEST COMPLETATI; MANCA LA RUN VERA NOTTURNA** | • Benchmark HPC eseguito: **`namd3 +p8` = `12.68 ns/day`** (**+40% più veloce** di `+p16` grazie all'assenza di contesa su AVX-512).<br>• Script d'analisi [`03_analyze_wt_baseline.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/03_analyze_wt_baseline.tcl) (che unisce `analysisMD.tcl` del prof + $Q$-value di Bekker) e [`plot_wt_baseline.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/plot_wt_baseline.py) validati e funzionanti.<br>• **ATTENZIONE IMPORTANTE:** I file `.dat`, i grafici `fig5`/`fig5b` e `unstable_residues_report.txt` presenti ora in `02_wt_baseline/` derivano **solo dal micro-test di 10 ps (`smoke_400K.dcd`, 11 frame)**! La simulazione vera (`wt_300K.conf` e `wt_400K.conf` tramite `run_overnight_wt.sh`) **non è ancora stata lanciata**. |
| **`03_mutations/`** | **Fase 3** (`ROAD3.md`): Screening Mutazioni & Scelta Mutanti | 🔄 **PRE-SCREENING STATICO FATTO (~70%); MANCA SCELTA FINALE + `BuildModel`** | • Allineati i 4 Ankyrin Repeats al consenso di Mosavi et al. (`04_ankyrin_consensus_and_targets.py`, grafico `fig6` generato).<br>• Eseguito `RepairPDB` (`2MYO_Repair.pdb`) e `PositionScan` su **43 posizioni candidate** (`PS_2MYO_Repair_scanning_output.txt`, `1.075` righe totali valutate), generando `foldx_top_mutations.txt` e `fig7`.<br>• Presenti ~1.000 file `.pdb` intermedi di `PositionScan` da pulire.<br>• Mancano: incrocio con i risultati della MD vera di stanotte, eventuali risultati dai server web (**DynaMut2 / FireProt 2.0 / ThermoMPNN-D**), scelta finale delle 2 mutazioni singole + 1 doppia, e `FoldX BuildModel`. |
| **`04_md_simulations/`** | **Fase 4** (`ROAD4.md`): MD Comparative (`WT` vs `Mut1` vs `Mut2` vs `Double`) | ⏳ **DA FARE DOPO LA FASE 3** | Cartella vuota: ospiterà i sistemi solvatati dei 3 mutanti e le 3 repliche *multi-seed* a `400 K`. |
| **`05_analysis/` & `06_presentation/`** | **Fase 5** (`ROAD5.md`): Analisi Statistica Finale & Slide (16 Dicembre) | ⏳ **PARZIALE (GIÀ PRESENTI I PRIMI PDF DI PROVA)** | `06_presentation/` contiene `fig5`/`fig5b` (da aggiornare dopo la run notturna vera), `fig6_ankyrin_consensus_alignment.pdf` (definitiva) e `fig7_foldx_position_scan_heatmap.pdf`. |

---

## 2. 🧠 La Logica Scientifica del Progetto (Il Filo Conduttore)

Prima di vedere i comandi, fissiamo **perché** ogni passo viene fatto esattamente in questo ordine (seguendo [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) Slide 15 e *Bekker et al., 2019*):

1. **Perché simuliamo il Wild-Type a `300 K` e `400 K` *prima* di scegliere definitivamente i mutanti?**  
   FoldX è un calcolatore statico (tiene lo scheletro proteico quasi rigido): ci dice *quanto cambia l'energia libera ($\Delta\Delta G$)* se sostituiamo un amminoacido, ma **non sa quali zone della proteina cedono davvero quando l'acqua scalda**.  
   Seguendo *Bekker et al. (2019)*:
   - A **`300 K` ($27^\circ\text{C}$, controllo):** la proteina nativa `2MYO` deve restare compatta e stabile ($Q \approx 0.90$, basso RMSD) $\rightarrow$ dimostra che il nostro setup `CHARMM36m` è sano.
   - A **`400 K` ($127^\circ\text{C}$, stress termico):** l'agitazione accelera i cedimenti locali senza distruggere il fold. Misurando **RMSF** e **frazione di contatti nativi per residuo ($\langle Q_i \rangle_{400\text{K}}$)** identifichiamo i **"bulloni deboli"** (residui con $\langle Q_i \rangle < 0.60$, o forte calo nei repeat/loop).
2. **Perché incrociamo MD (`<Q_i>_400K`) + FoldX (`PositionScan`) + Consenso Ankyrin + AI/NMA (ThermoMPNN/FireProt/DynaMut2)?**  
   Per trasformare il progetto in un lavoro da ricerca (*Orthogonal Consensus Funnel*): scegliamo le **2 mutazioni singole (`Mut1`, `Mut2`)** e il **Doppio Mutante (`Mut1 + Mut2`)** solo dove **tutti i metodi convergono**:
   - **Fuori dalla Functional Mask** (i 24 residui di `3AAA` restano intoccabili $\rightarrow$ funzione biologica preservata).
   - **Zona vulnerabile nella MD a `400 K`** ($\langle Q_i \rangle$ basso o picco di RMSF $\rightarrow$ ha senso strutturale rinforzarla).
   - **Energia favorevole in FoldX** ($\Delta\Delta G_{\text{FoldX}} < -1.0\text{ kcal/mol}$) **e meccanismo chimico chiaro** (es. miglior impaccamento idrofobico $\Delta E_{\text{VdW+SolvH}}$, rilascio di clash/tensione di Ramachandran, nuovi H-bond).
   - **Conferma ortogonale** (Consenso Ankyrin di *Mosavi et al., 2002* e/o predittori *ThermoMPNN-D / FireProt 2.0 / DynaMut2*).
3. **Perché poi simuliamo `WT`, `Mut1`, `Mut2` e `Double` con più repliche (`seed` diversi) a `400 K`?**  
   Una singola traiettoria MD è stocastica (dipende dalle velocità iniziali). Ripetendo le simulazioni a `400 K` con **3 `seed` indipendenti** per tutte e 4 le varianti con identico protocollo (*Matched Comparison*, Slide 21), calcoliamo **Media $\pm$ Deviazione Standard ($\mu \pm \sigma$)** di $Q$ (globale, per classi chimiche e per singolo `ANK1–4`) e dimostriamo statisticamente che i nostri mutanti trattengono più contatti nativi del Wild-Type!

---

## 3. 📝 Scaletta Cronologica Passo-Passo (`STEP 1` $\rightarrow$ `STEP 10`)

### 🟢 FASE 1 — Analisi Strutturale & Setup Scatola d'Acqua (✅ COMPLETATA)

- **STEP 1 — Ispezione Biofisica & Functional Mask (`2MYO` vs `3AAA`) [✅ Fatto]**
  - *Cosa:* Analizzati i record PDB, Cisteine, Istidine e l'interfaccia $< 4.5\text{ Å}$ con Actin-Capping Protein (`3AAA`).
  - *Perché:* Garantire 0 gap, tautomeri corretti (`HIS 40 -> HSD`, `HIS 67 -> HSE`) e proteggere i 24 residui funzionali.
  - *Script:* [`scripts/01_inspect_2myo.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/01_inspect_2myo.tcl) + foto VMD [`viz_01_structure_and_mask.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_01_structure_and_mask.tcl).

- **STEP 2 — Topologia CHARMM36m, Solvatazione `TIP3P` e Ionizzazione `0.15 M NaCl` [✅ Fatto]**
  - *Cosa:* Rimozione idrogeni NMR, `psfgen`, box d'acqua `12 Å` e neutralizzazione salina (`21.461` atomi).
  - *Perché:* Creare il sistema fisico periodico su cui girerà la Dinamica Molecolare.
  - *Output:* `01_structures/2myo_wb_ionized.psf`, `.pdb`, `pbc_box_params.txt` + foto VMD [`viz_02_solvated_box.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_02_solvated_box.tcl).

---

### 🟡 FASE 2 — Benchmark HPC & Simulazione Baseline del Wild-Type

- **STEP 3 — Benchmark HPC (`+p8` vs `+p16`) & Validazione Script d'Analisi [✅ Fatto]**
  - *Cosa:* Misurata la velocità reale NAMD 3 (`12.68 ns/day` su 8 core fisici `+p8`) e validati [`scripts/03_analyze_wt_baseline.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/03_analyze_wt_baseline.tcl) + [`scripts/plot_wt_baseline.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/plot_wt_baseline.py) sullo smoke-test da 10 ps.
  - *Perché:* Evitare sprechi di tempo CPU (speedup +40% usando `+p8` invece di `+p16`) ed essere certi che l'analisi automatica non crashi dopo 7 ore di calcolo.

- **STEP 4 — Simulazione MD Vera del Wild-Type (`300 K` Controllo + `400 K` Stress Termico) [🌙 PRIORITÀ #1: DA LANCIARE STANOTTE]**
  - *Cosa:* Eseguire la run di produzione notturna per `WT 300 K` (`1.5 ns` = `750.000` step) e `WT 400 K` (`2.5 ns` = `1.250.000` step, ~7.5 ore totali; oppure estesa a `5 ns` se lasciata più ore), seguita dall'analisi automatica VMD + PDF.
  - *Perché:* Ottenere la **vera dinamica di equilibrio e cedimento termico** del Wild-Type (`wt_300K.dcd` e `wt_400K.dcd`), sostituendo i dati provvisori dello smoke-test da 10 ps con i veri valori di **RMSD**, **RMSF**, **H-Bonds**, **SASA**, **$R_g$**, **% $\alpha$-elica**, **$Q(t)$ per Repeat (`ANK1–4`)** e **$\langle Q_i \rangle_{400\text{K}}$ per singolo residuo**.
  - *Comando (dalla cartella `BML/`):*
    ```bash
    nohup bash 02_wt_baseline/run_overnight_wt.sh 8 > 02_wt_baseline/overnight_master.log 2>&1 &
    ```
  - *Output atteso al risveglio:* `02_wt_baseline/wt_300K.dcd`, `wt_400K.dcd`, `unstable_residues_report.txt` (definitivo), `wt_400K_q_beta.pdb`, `06_presentation/fig5_wt_baseline_plots.pdf`, `fig5b_wt_analysisMD_observables.pdf`.

- **STEP 5 — Ispezione Risultati Reali WT & Foto 3D Instabilità [☀️ DOMANI MATTINA]**
  - *Cosa:* Controllare `02_wt_baseline/unstable_residues_report.txt`, verificare l'equilibrio nei PDF `fig5`/`fig5b`, e lanciare la scena 3D di VMD per fotografare i punti deboli a `400 K`.
  - *Perché:* Identificare visivamente e quantitativamente quali Ankyrin Repeats e quali amminoacidi perdono contatti a `400 K`.
  - *Comandi:*
    ```bash
    cat 02_wt_baseline/unstable_residues_report.txt
    vmd -e scripts/viz_03_wt_instability.tcl
    ```

---

### 🟠 FASE 3 — Screening Multi-Tool & Costruzione 3D dei Mutanti

- **STEP 6 — Screening Statico & Matrice di Consenso Multi-Tool [🔄 70% GIÀ FATTO; INCROCIO FINALE DOMANI]**
  - *Regole e Motivazione delle Scelte:* Vedi il documento dedicato **[`MUTATION_SELECTION.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/MUTATION_SELECTION.md)** (spiega le 6 Regole oggettive di inclusione/esclusione e come **MD 400K + FireProt/Rosetta + ThermoMPNN-D + Consenso Ankyrin + FoldX + DynaMut2** convergono senza scelte arbitrarie).
  - *Cosa già fatta:*
    - Consenso Ankyrin (`scripts/04_ankyrin_consensus_and_targets.py` $\rightarrow$ `06_presentation/fig6_ankyrin_consensus_alignment.pdf`).
    - `FoldX RepairPDB` (`2MYO_Repair.pdb`) e `PositionScan` su **43 posizioni** (`PS_2MYO_Repair_scanning_output.txt` + `foldx_top_mutations.txt` + `fig7`).
    - Preparato [`scripts/06_build_consensus_matrix.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/06_build_consensus_matrix.py) e [`03_mutations/dynamut2_mutations_list.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/dynamut2_mutations_list.txt).
  - *Cosa fare domani dopo lo STEP 5:* Rilanciare `python3 scripts/06_build_consensus_matrix.py` (che leggerà i veri `<Q_i>_400K` della MD + FoldX + eventuali CSV da FireProt/DynaMut2/ThermoMPNN-D) e selezionare ufficialmente:
    1. **`Mut1`** (prima mutazione singola stabilizzante su zona debole MD + consenso)
    2. **`Mut2`** (seconda mutazione singola in sito/repeat distinto o complementare)
    3. **`Double`** (`Mut1 + Mut2`, combinazione doppia per testare additività vs epistasi)

- **STEP 7 — Modellazione 3D dei 3 Mutanti con FoldX `BuildModel` & Foto Zoom [☀️ DOMANI MATTINA — 2 minuti]**
  - *Cosa:* Scrivere `03_mutations/individual_list.txt` con le 3 varianti (`Mut1;`, `Mut2;`, `Mut1,Mut2;`) ed eseguire `foldx --command=BuildModel` (`numberOfRuns=3..5`).
  - *Perché:* Generare le coordinate atomiche 3D dei 3 mutanti (`2MYO_Repair_1.pdb`, `_2.pdb`, `_3.pdb`) ed estrarre la scomposizione energetica ufficiale `Dif_2MYO_Repair.fxout` per spiegare chimicamente ogni mutazione nelle slide. Creare poi lo script VMD `viz_04_mutations_zoom.tcl` per fotografare il prima/dopo della tasca mutata.

---

### 🔵 FASE 4 — Simulazioni MD Comparative ("Matched Comparison")

- **STEP 8 — Preparazione Identica dei Sistemi Solvatati per i 3 Mutanti [☀️ DOMANI — 2 minuti]**
  - *Cosa:* Creare uno script TCL parametrico in `scripts/` che applica ai 3 PDB mutati usciti da FoldX esattamente la stessa pipeline `CHARMM36m` del Wild-Type (`psfgen` con `HIS 40 -> HSD`, `HIS 67 -> HSE`, `solvate -t 12`, `autoionize -sc 0.15`, stesse dimensioni box).
  - *Perché:* Obbligo assoluto di rigore scientifico ("Matched Comparison", `Laboratory_Project_2026.pdf` Slide 3 & 21): le uniche differenze tra le simulazioni devono essere le mutazioni (e non il numero d'atomi o i parametri del box).

- **STEP 9 — Campagna MD Multi-Seed a `400 K` (`WT`, `Mut1`, `Mut2`, `Double`) [🔥 NOTTI SUCCESSIVE]**
  - *Cosa:* Eseguire **almeno 3 repliche indipendenti** (variando `seed` delle velocità di Maxwell-Boltzmann) a `400 K` per ciascuno dei **4 sistemi**:
    - `WT` a `400 K`: Replica 1 *(già fatta stanotte nello STEP 4!)* + Replica 2 + Replica 3
    - `Mut1` a `400 K`: Replica 1 + Replica 2 + Replica 3
    - `Mut2` a `400 K`: Replica 1 + Replica 2 + Replica 3
    - `Double` a `400 K`: Replica 1 + Replica 2 + Replica 3
  - *Suggerimento pratico:* Sono 11 run rimanenti a `400 K`: possiamo impostare un batch script per la notte o **dividerle a metà tra il tuo PC (`+p8`, Linux) e il PC di Alessia (Windows)** dimezzando i tempi!

---

### 🟣 FASE 5 — Analisi Statistica Multi-Replica, Epistasi & Presentazione

- **STEP 10 — Analisi Comparativa Finale, Grafici `.pdf` & Slide (16 Dicembre) [📊 FASE FINALE]**
  - *Cosa:*
    1. Calcolare su tutte le repliche (finestra finale equilibrata, **Media $\pm$ SD** e **$p$-value / $t$-test**):
       - **$Q$-value** (`all-all`, `hydrophilic-all`, `hydrophobic-small` — *Bekker et al. Table II*) e **$Q$ per singolo Repeat (`ANK1`..`ANK4`)**.
       - Conteggio dei **residui instabili residui ($\langle Q_i \rangle < 0.60$)** per variante.
       - **$\Delta\text{RMSF}$**, **H-bond occupancy**, **SASA**, **$R_g$** e **% di $\alpha$-elica** (`analysisMD.tcl`).
       - **Ciclo Epistatico del Doppio Mutante** ($\Delta Q_{1+2}$ vs $\Delta Q_1 + \Delta Q_2$) e mappa differenza contatti $\Delta C_{ij}$.
    2. Esportare tutti i grafici in **`.pdf` vettoriale** in `06_presentation/` e montare le slide (max 10 min) con la divisione dei contributi di Oleg e Alessia.

---

## 4. ⚡ Cosa Possiamo Fare ADESSO (Oggi Pomeriggio, 15 Minuti) per Anticipare e Migliorare il Lavoro

Prima di lanciare la simulazione notturna nello **STEP 4**, puoi fare queste **4 azioni rapide e super logiche** che mettono in ordine la repo e fanno lavorare gratis anche i cluster accademici remoti durante la notte:

### Azione 1 — Pulire i ~1.000 PDB temporanei da `03_mutations/` (5 secondi)
Durante `PositionScan`, FoldX ha creato circa ~1.000 file intermedi come `ALA24_2MYO_Repair.pdb`, `VAL110_2MYO_Repair.pdb` ecc. che sporcano la vista della cartella senza servirci più (i risultati energetici sono già salvati in `PS_2MYO_Repair_scanning_output.txt` e `energies_*_2MYO_Repair.txt`).  
Per fare pulizia **senza toccare `2MYO.pdb`, `2MYO_Repair.pdb` né i report**:
```bash
# Controlla prima quanti sono (~1000 file) e poi rimuovi solo i PDB temporanei dei singoli residui:
ls 03_mutations/[A-Z]*[0-9]*_2MYO_Repair.pdb | wc -l
rm -f 03_mutations/[A-Z]*[0-9]*_2MYO_Repair.pdb
```

### Azione 2 — Generare subito la Matrice di Scomposizione Energetica FoldX (2 secondi)
Anche se i `<Q_i>_400K` definitivi arriveranno domani mattina, puoi già far girare [`scripts/06_build_consensus_matrix.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/06_build_consensus_matrix.py) per vedere la **scomposizione fisica** ($\Delta E_{\text{VdW+SolvH}}$, $\Delta E_{\text{Clash}}$, $\Delta E_{\text{HBond}}$) delle tue migliori mutazioni FoldX ed evitare falsi positivi (es. Cisteine libere o Proline dentro $\alpha$-eliche):
```bash
python3 scripts/06_build_consensus_matrix.py
cat 03_mutations/consensus_master_table.txt
```
Questo genera anche il nuovo grafico vettoriale **`06_presentation/fig8_consensus_energy_breakdown.pdf`**!

### Azione 3 — Sottomettere i Job ai Server Web Gratuiti (Compito di Alessia — vedi [`ROAD_PER_ALE/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROAD_PER_ALE))
Mentre il tuo PC fa la MD stanotte, Alessia può far lavorare in parallelo i server remoti seguendo la guida passo-passo **[`ROAD_PER_ALE/01_GUIDA_OPERATIVA_ALESSIA.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROAD_PER_ALE/01_GUIDA_OPERATIVA_ALESSIA.md)** (e vedendo in **[`ROAD_PER_ALE/02_LAVORO_OLEG_E_UNIONE_DATI.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROAD_PER_ALE/02_LAVORO_OLEG_E_UNIONE_DATI.md)** come uniremo i dati domani mattina) usando **`03_mutations/2MYO_Repair.pdb`**:
1. **FireProt 2.0** ([https://loschmidt.chemi.muni.cz/fireprot/](https://loschmidt.chemi.muni.cz/fireprot/)): carica `03_mutations/2MYO_Repair.pdb` (Chain `A`), inserisci la mail `@mail.polimi.it` e invia (farà girare **Rosetta + FoldX + Filogenesi MSA** sui loro server!).
2. **DynaMut2** ([https://biosig.lab.uq.edu.au/dynamut2/](https://biosig.lab.uq.edu.au/dynamut2/) $\rightarrow$ **`Run > Single Mutation` $\rightarrow$ `Mutation List`**): carica `03_mutations/2MYO_Repair.pdb` e la lista già pronta [`03_mutations/dynamut2_mutations_list.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/dynamut2_mutations_list.txt) per avere i $\Delta\Delta G$ basati su *Normal Mode Analysis*.
3. **ThermoMPNN-D Colab** ([Apri Colab ThermoMPNN-D](https://colab.research.google.com/github/Kuhlman-Lab/ThermoMPNN-D/blob/main/ThermoMPNN-D.ipynb)): carica `03_mutations/2MYO_Repair.pdb` su GPU Google gratuita per ottenere le predizioni AI su tutte le mutazioni singole ed epistasi doppie.
*(Appena Alessia mette i `.csv` in `03_mutations/`, basterà rilanciare `python3 scripts/06_build_consensus_matrix.py`!).*

### Azione 4 — Impostare la Durata Ideale e Lanciare la MD Notturna del Wild-Type (Stasera)
Attualmente:
- [`02_wt_baseline/wt_300K.conf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/wt_300K.conf) è impostato a **`1.5 ns`** (`750.000` step $\rightarrow$ **~2.8 ore** su `+p8`).
- [`02_wt_baseline/wt_400K.conf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/wt_400K.conf) è impostato a **`2.5 ns`** (`1.250.000` step $\rightarrow$ **~4.7 ore** su `+p8`).
- **Totale attuale:** **~7.5 ore** (perfetto per una singola notte!).  
*(Consiglio opzionale per migliorare ancora di più la statistica a `400 K` senza rischiare il tempo: se stanotte il PC può girare ~10.5–12.5 ore, puoi portare `run_steps` in `wt_400K.conf` da `1250000` (`2.5 ns`) a `2000000` (`4.0 ns`, ~10.4 ore totali) o `2500000` (`5.0 ns`, ~12.3 ore totali); altrimenti lascia i valori attuali per trovarlo sicuramente finito dopo ~7.5 ore!)*

Quando sei pronto per andare a dormire, lancia dalla cartella `BML/`:
```bash
nohup bash 02_wt_baseline/run_overnight_wt.sh 8 > 02_wt_baseline/overnight_master.log 2>&1 &
```
Per sbirciare l'avanzamento in qualsiasi momento senza fermarlo:
```bash
tail -f 02_wt_baseline/overnight_master.log
```
*(Premi `Ctrl+C` per uscire dal `tail`: NAMD continuerà a girare in background).*
