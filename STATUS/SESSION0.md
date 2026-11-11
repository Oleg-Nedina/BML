# 📅 SESSION0.md — Diario della Sessione 0 (Setup, Fase 1 Completa e Preparazione Batch Fase 2)

**Data:** 11 Novembre 2026  
**Partecipanti:** Oleg Nedina & Alessia Rigoni (Team #10)  
**Stato Avanzamento:** ✅ **Punto 0 Completato** | ✅ **Fase 1 (`ROAD1.md`) Completata al 100%** | 🔄 **Fase 2 (`ROAD2.md`) Validata e Pronta per il Run Notturno** | 🚀 **Fase 3 (`ROAD3.md`) Impostata e Avviata in Anticipo!**

---

## 🎯 1. Obiettivo della Sessione di Oggi

1. Configurare e verificare l'intero ambiente software (**NAMD 3.0.3**, **VMD 2.0.0** con `VMDextensions.tcl`, **FoldX 5** e **Python 3**) e strutturare la repository `BML` con regole rigorose da lavoro di ricerca ([`AGENTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/AGENTS.md)).
2. Svolgere interamente la **Fase 1 ([`ROADMAP/ROAD1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/ROAD1.md))**: ispezione strutturale e biofisica di `2MYO.pdb`, verifica delle Cisteine e delle Istidine, estrazione oggettiva della *Functional Mask* dal complesso cristallografico `3AAA.pdb`, e costruzione del dominio di simulazione solvatato e ionizzato in `CHARMM36m` (`01_structures/2myo_wb_ionized.psf` / `.pdb`).
3. Impostare, integrare con `analysisMD.tcl` del Prof. Gautieri e validare tramite *smoke test* l'intera pipeline della **Fase 2 ([`ROADMAP/ROAD2.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/ROAD2.md))**, lasciando pronto lo script batch per l'esecuzione notturna delle simulazioni MD del Wild-Type a **$300\text{ K}$** e **$400\text{ K}$**.

---

## 🛠️ 2. Cosa Abbiamo Fatto Oggi (Riepilogo Passo-Passo)

### A. Setup dell'Ambiente e Regole di Progetto (Punto 0)
- **Installazione di FoldX 5:** Estratto `foldx5_Linux_0.zip` in `/usr/local/foldx/foldx_20261231` e creato il symbolic link `/usr/local/bin/foldx`.
- **Verifica Toolchain:** Verificati `namd3`, `vmd` (con `prepareNativeContacts` da `VMDextensions.tcl`), `foldx` e le librerie Python (`numpy`, `scipy`, `pandas`, `matplotlib`).
- **Portabilità Linux/Windows:** Scaricati `2MYO.pdb` e `3AAA_complex_ref.pdb` in `01_structures/raw/` e copiati i file del force field `CHARMM36m` (`top_all36_prot.rtf`, `par_all36m_prot.prm`, `toppar_water_ions.str`) in `common/charmm36/` e `VMDextensions.tcl` in `scripts/`, così tutti gli script usano percorsi relativi identici su Linux e Windows.
- **Definizione Regole ([`AGENTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/AGENTS.md)):** Fissate le regole del progetto: *Zero Ipotesi*, *Citazione obbligatoria delle fonti*, *Note Foto VMD per le slide*, *Grafici obbligatoriamente in formato vettoriale `.pdf`* e tracciamento continuo in `ROADMAP/` e `STATUS/`.

### B. Completamento della Fase 1 (`ROADMAP/ROAD1.md` — Studio di `2MYO` e Costruzione del Dominio)
1. **Ispezione Chimico-Fisica ([`scripts/01_inspect_2myo.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/01_inspect_2myo.tcl)):**
   - Verificata la completezza della sequenza: **118 amminoacidi continui** (`902` atomi pesanti + `913` idrogeni NMR = `1815` atomi).
   - Verificate le **3 Cisteine (`CYS 2, 45, 83`)**: distanze $S_\gamma\text{--}S_\gamma$ pari a `21.88 Å`, `27.01 Å` e `9.43 Å` $\rightarrow$ **zero ponti disolfuro** (tutte ridotte).
   - Determinato oggettivamente lo stato di protonazione delle **4 Istidine**:
     - **`HIS 40` $\rightarrow$ `HSD`** (il protone su `ND1` dona un H-bond inter-repeat a `SER 74 : OG` a $2.44\text{ Å}$).
     - **`HIS 67` $\rightarrow$ `HSE`** (`ND1` accetta un H-bond da `LYS 37 : NZ` a $3.24\text{ Å}$, mentre `NE2-HE2` dona un H-bond a `GLY 35 : O` a $2.40\text{ Å}$).
     - **`HIS 68` e `HIS 80` $\rightarrow$ `HSD`** (esposte/apolari, mantengono l'assegnazione NMR).
   - Misurato l'**RMSD tra `2MYO` (Ratto, NMR) e `3AAA` (`Chain C`, Uomo, X-Ray $2.20\text{ Å}$)**: **`2.365 Å`**, ed estratti i **24 residui della Functional Mask** a distanza $< 4.5\text{ Å}$ dall'Actin-Capping Protein (`resid 4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100`).
2. **Costruzione Topologia `CHARMM36m`, Solvatazione `TIP3P` e Ionizzazione `0.15 M NaCl` ([`scripts/02_build_wt_charmm36m.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/02_build_wt_charmm36m.tcl)):**
   - Generati `01_structures/2myo_wt.psf` e `.pdb` (carica netta della proteina: **$-5.00\text{ e}$**, `116` termini CMAP cross-term attivi).
   - Solvatata in box `TIP3P` con padding di $12.0\text{ Å}$ e ionizzata a $0.15\text{ M NaCl}$ (`01_structures/2myo_wb_ionized.psf` e `.pdb`):
     - **Atomi totali:** `21.461` (`6.535` molecole d'acqua, `23` ioni $\text{Na}^+$, `18` ioni $\text{Cl}^-$, carica netta finale = **`0.0000 e`**).
     - **Dimensioni cella PBC:** `71.840 × 58.972 × 54.840 Å`, centro `(7.735, -8.087, 86.003)` (salvate in `01_structures/pbc_box_params.txt`).
3. **Script Grafici 3D per le Slide della Fase 1:**
   - Creati e testati [`scripts/viz_01_structure_and_mask.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_01_structure_and_mask.tcl) e [`scripts/viz_02_solvated_box.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_02_solvated_box.tcl).

### C. Benchmark HPC, Integrazione di `analysisMD.tcl` e Validazione della Fase 2 (`ROADMAP/ROAD2.md`)
1. **Benchmark HPC (`+p16` vs `+p8`):**
   - Eseguito il benchmark su `21.461` atomi scoprendo che **`namd3 +p8` (8 core fisici puri)** raggiunge **`12.68 ns/day`** (`0.0133 s/step`), battendo del **$+40\%$** l'esecuzione su 16 thread SMT (`+p16`, `9.01 ns/day`) grazie all'assenza di contesa sulle unità vettoriali AVX-512!
2. **Integrazione di `analysisMD.tcl` del Prof. Gautieri + $Q$-value di Bekker et al.:**
   - Creato [`scripts/03_analyze_wt_baseline.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/03_analyze_wt_baseline.tcl) che unisce tutte le metriche di `laboratorio/analysisMD.tcl` (**Allineamento sul backbone, RMSD, RMSF sui $C_\alpha$, Hydrogen Bonds + dettagli di occupancy, SASA, Raggio di Girazione $R_g$, % di $\alpha$-elica STRIDE**) con il calcolo del **$Q$-value** (`all-all`, `hydrophilic-all`, `hydrophobic-small`, per singolo Repeat `ANK1..4`, e $\langle Q_i \rangle$ per ciascuno dei 118 residui).
   - Creato [`scripts/plot_wt_baseline.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/plot_wt_baseline.py) che esporta i due pannelli a 4 quadranti in **formato vettoriale `.pdf`** (`06_presentation/fig5_wt_baseline_plots.pdf` e `fig5b_wt_analysisMD_observables.pdf`).
   - Creata la documentazione completa di tutti gli script in [`scripts/README_SCRIPTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/README_SCRIPTS.md).

### D. Preparazione Anticipata della Fase 3 ([`ROADMAP/ROAD3.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/ROAD3.md))
Per massimizzare il tempo prima di lanciare la simulazione MD notturna, abbiamo già strutturato la **Fase 3 (*Imbuto Decisionale Multi-Tool*)**:
1. **Funnel di Selezione delle Posizioni (`118` $\rightarrow$ `94` $\rightarrow$ `88` $\rightarrow$ `35` Hotspot $\rightarrow$ `2` Mutazioni Singole + `1` Doppia):**
   - **Da 118 a 94 residui:** Esclusi i **24 residui della Functional Mask** (`3AAA` $< 4.5\text{ Å}$).
   - **Da 94 a 88 residui interni:** Escluse le code terminali libere (`resid 1–3` e `116–118`).
   - **Pre-screening di 35 Hotspot Strutturali e di Consenso ([`scripts/04_ankyrin_consensus_and_targets.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/04_ankyrin_consensus_and_targets.py)):** Allineati i 4 Ankyrin Repeats al consenso universale di Mosavi et al. (PNAS 2002) ed estratti tutti i Ramachandran outliers dal `REMARK 500` di `2MYO.pdb`, generando il grafico vettoriale `06_presentation/fig6_ankyrin_consensus_alignment.pdf` e il file `03_mutations/position_scan_hotspots.cfg` (`35 posizioni × 20 aa = 700 mutazioni`).
   - **Configurazione FoldX (`03_mutations/repair.cfg` + `scripts/05_parse_foldx_screening.py`):** Pronti per eseguire `RepairPDB` e `PositionScan` già oggi pomeriggio, più invio in Cloud a **ThermoMPNN** (che valuta tutte le `118 × 20 = 2.360` mutazioni in 15 secondi) e **FireProt 2.0**, così domani mattina basterà incrociare la tabella energetica con i residui instabili ($\langle Q_i \rangle_{400\text{K}} < 0.60$) della MD notturna!

---

## 🌙 3. Cosa Fare Stanotte (Lancio del Batch Notturno della Fase 2)

Prima di andare a dormire, dalla cartella principale `/home/olly/UNI/BioMolEng2026/laboratorio/BML`:

1. **(Consigliato prima di spegnere lo schermo):** Rilancia i 2 comandi di validazione rapida (15 secondi) per verificare che i due PDF di test vengano creati:
   ```bash
   vmd -dispdev text -e scripts/03_analyze_wt_baseline.tcl
   python3 scripts/plot_wt_baseline.py
   ```
2. **Quando vuoi far partire la simulazione vera per la notte (~7.5 ore totali su `+p8`):**
   ```bash
   nohup bash 02_wt_baseline/run_overnight_wt.sh 8 > 02_wt_baseline/overnight_master.log 2>&1 &
   ```
   *(Con `nohup ... &` il processo continua a girare in background anche se chiudi il terminale, e puoi controllare a che punto è in qualsiasi momento con `tail -f 02_wt_baseline/overnight_master.log` o `tail -f 02_wt_baseline/wt_300K.log`).*

Al tuo risveglio domani mattina, lo script avrà già eseguito da solo:
- La run **WT a 300 K** (`1.5 ns`, file `02_wt_baseline/wt_300K.dcd`),
- La run **WT a 400 K** (`2.5 ns`, file `02_wt_baseline/wt_400K.dcd`),
- Tutta l'analisi VMD (`02_wt_baseline/unstable_residues_report.txt` e tutti i `.dat`),
- I due grafici vettoriali PDF definitivi (`06_presentation/fig5_wt_baseline_plots.pdf` e `fig5b_wt_analysisMD_observables.pdf`)!

---

## ⏭️ 4. Programma per la Prossima Sessione (`SESSION1.md`)

Nella prossima sessione (appena il run notturno avrà finito):
1. **Analisi dei Risultati Reali della Fase 2 (`02_wt_baseline/unstable_residues_report.txt` e Grafici PDF):**
   - Leggeremo insieme quali dei 4 Ankyrin Repeats (`ANK1..4`) ha perso più contatti a $400\text{ K}$ e quali sono i "bulloni instabili" ($\langle Q_i \rangle < 0.60$) emersi dalla simulazione vera.
   - Scatteremo la **Foto 3D dell'instabilità a 400 K** con `vmd -e scripts/viz_03_wt_instability.tcl`.
2. **Generazione ed Esecuzione di `ROADMAP/ROAD3.md` (Fase 3 — Multi-Tool Consensus Screening & Modellazione dei Mutanti):**
   - Eseguiremo `RepairPDB` e `PositionScan` di **FoldX 5** sui residui instabili emersi dalla MD e sulle posizioni che deviano dal consenso degli Ankyrin Repeats (`LEU 73`, `LYS 94`, ecc.).
   - Incroceremo i $\Delta\Delta G$ di **FoldX** con **FireProt 2.0 / Rosetta / ThermoMPNN / DynaMut2** per selezionare i **2 migliori mutanti singoli (`Mut1`, `Mut2`)** e il **Doppio Mutante (`Mut1 + Mut2`)**.
   - Costruiremo i modelli 3D dei mutanti con `BuildModel` di FoldX e prepareremo le scatole d'acqua per le simulazioni comparative della Fase 4!
