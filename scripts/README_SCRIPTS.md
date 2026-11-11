# 🛠️ Guida Completa agli Script del Progetto (`scripts/README_SCRIPTS.md`)

Questo documento spiega in modo semplice e operativo **come lanciare ogni script** presente nella repository `BML`, **cosa fa ciascuno script** e **quale output (testuale, numerico e visivo 3D/2D) dovete aspettarvi** per verificare che tutto stia funzionando correttamente.

> [!IMPORTANT]
> **Regola base di esecuzione:** Tutti i comandi vanno lanciati stando nella cartella principale della repository (`/home/olly/UNI/BioMolEng2026/laboratorio/BML` su Linux o `...\BML` su Windows), perché gli script usano percorsi relativi a partire da `BML/`.

---

## 📋 Indice Rapido degli Script

| # | Script / File | Fase | Tipo | Comando per lanciarlo |
| :---: | :--- | :---: | :---: | :--- |
| **0** | `scripts/VMDextensions.tcl` | Tutte | Libreria TCL | *(Caricata in automatico dagli altri script)* |
| **1** | `scripts/01_inspect_2myo.tcl` | Fase 1 | Analisi VMD | `vmd -dispdev text -e scripts/01_inspect_2myo.tcl` |
| **2** | `scripts/viz_01_structure_and_mask.tcl` | Fase 1 | **Grafica 3D VMD** | `vmd -e scripts/viz_01_structure_and_mask.tcl` |
| **3** | `scripts/02_build_wt_charmm36m.tcl` | Fase 1 | Setup Topologia | `vmd -dispdev text -e scripts/02_build_wt_charmm36m.tcl` |
| **4** | `scripts/viz_02_solvated_box.tcl` | Fase 1 | **Grafica 3D VMD** | `vmd -e scripts/viz_02_solvated_box.tcl` |
| **5** | `02_wt_baseline/smoke_test_400K.conf` | Fase 2 | Test NAMD3 | `namd3 +p8 02_wt_baseline/smoke_test_400K.conf > 02_wt_baseline/smoke_400K.log` |
| **6** | `scripts/03_analyze_wt_baseline.tcl` | Fase 2 | Analisi MD + `analysisMD.tcl` | `vmd -dispdev text -e scripts/03_analyze_wt_baseline.tcl` |
| **7** | `scripts/plot_wt_baseline.py` | Fase 2 | **Grafici 2D Python** | `python3 scripts/plot_wt_baseline.py` |
| **8** | `scripts/viz_03_wt_instability.tcl` | Fase 2 | **Grafica 3D VMD** | `vmd -e scripts/viz_03_wt_instability.tcl` |
| **9** | `02_wt_baseline/run_overnight_wt.sh` | Fase 2 | **Batch Notturno** | `bash 02_wt_baseline/run_overnight_wt.sh` |
| **10** | `scripts/04_ankyrin_consensus_and_targets.py` | Fase 3 | Consenso Ankyrin + Config FoldX | `python3 scripts/04_ankyrin_consensus_and_targets.py` |
| **11** | `03_mutations/repair.cfg` | Fase 3 | FoldX `RepairPDB` | `foldx -f 03_mutations/repair.cfg > 03_mutations/repair.log` |
| **12** | `03_mutations/position_scan_hotspots.cfg` | Fase 3 | FoldX `PositionScan` | `foldx -f 03_mutations/position_scan_hotspots.cfg > 03_mutations/position_scan.log` |
| **13** | `scripts/05_parse_foldx_screening.py` | Fase 3 | **Classifica + Grafico PDF FoldX** | `python3 scripts/05_parse_foldx_screening.py` |

---

## 🔍 Dettaglio di Ogni Script e Output Atteso (Numerico + Visivo)

---

### 1. `scripts/01_inspect_2myo.tcl` — Ispezione Chimico-Fisica di `2MYO` e `3AAA`
- **Come si lancia:**
  ```bash
  vmd -dispdev text -e scripts/01_inspect_2myo.tcl
  ```
- **Cosa fa:** Legge `01_structures/raw/2MYO.pdb` e `3AAA_complex_ref.pdb` senza interfaccia grafica e misura la sequenza FASTA, le distanze tra le 3 Cisteine, i vicini polari delle 4 Istidine e i residui a contatto ($< 4.5\text{ Å}$) con l'Actin-Capping Protein.
- **Cosa aspettarsi come output (a schermo):**
  - `Numero totale atomi nel PDB: 1815` (`902` pesanti + `913` idrogeni), `118` residui continui senza `/`.
  - Distanze $S_\gamma\text{--}S_\gamma$ delle Cisteine: `21.88 Å`, `27.01 Å`, `9.43 Å` (nessun ponte disolfuro).
  - Contatti delle Istidine: `HIS 40` (`ND1` a `2.44 Å` da `SER 74 OG` $\rightarrow$ **`HSD`**) e `HIS 67` (`ND1` a `3.24 Å` da `LYS 37 NZ` e `NE2` a `2.40 Å` da `GLY 35 O` $\rightarrow$ **`HSE`**).
  - `RMSD tra 2MYO e 3AAA`: **`2.365 Å`** e lista dei **24 residui della Functional Mask** (`resid 4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100`).

---

### 2. `scripts/viz_01_structure_and_mask.tcl` — 📸 Visualizzazione 3D dei 4 Repeat e della Functional Mask (Slide 2 & 3)
- **Come si lancia:**
  ```bash
  vmd -e scripts/viz_01_structure_and_mask.tcl
  ```
- **Cosa aspettarsi come output VISIVO (Finestra 3D di VMD):**
  - Si apre la finestra grafica di VMD con **sfondo bianco pulito** e senza assi cartesiani.
  - Vedrete una proteina a forma di "molla curva" (solenoide $\alpha$-elica) divisa in **4 colori**:
    - **Blu** = `ANK1` (residui 1–34, N-terminale)
    - **Ciano** = `ANK2` (residui 35–67)
    - **Arancione** = `ANK3` (residui 68–100)
    - **Verde** = `ANK4` (residui 101–118, C-terminale)
  - Sulla faccia concava dei loop vedrete **3 bastoncini Rossi** che sporgono: sono **`TRP 8`, `ASP 44`, `GLU 78`** (i 3 ganci funzionali).
  - Al centro vedrete i bastoncini colorati per atomo delle istidine **`HIS 40`** e **`HIS 67`** e dei loro partner di legame a idrogeno.
  - **Come fare le foto:** Aprendo `Graphics -> Representations`, potete fare doppio click sulla riga `Surf` per accendere una "nuvola" rossa semitrasparente sopra i 24 amminoacidi della zona funzionale vietata, e salvare l'immagine da `File -> Render -> Snapshot`!

---

### 3. `scripts/02_build_wt_charmm36m.tcl` — Costruzione Topologia `CHARMM36m`, Scatola d'Acqua e Ioni
- **Come si lancia:**
  ```bash
  vmd -dispdev text -e scripts/02_build_wt_charmm36m.tcl
  ```
- **Cosa fa:** Toglie i vecchi idrogeni NMR, imposta i tautomeri corretti delle 4 istidine, ricostruisce la topologia `CHARMM36m` (`psfgen`), immerge la proteina in un cubo d'acqua `TIP3P` con margine $12\text{ Å}$ (`solvate`) e aggiunge ioni $\text{Na}^+$ e $\text{Cl}^-$ a $0.15\text{ M}$ per portare la carica totale a zero (`autoionize`).
- **Cosa aspettarsi come output (File generati + Testo):**
  - **File creati in `01_structures/`:**
    - `2myo_wt.psf` e `2myo_wt.pdb` (solo proteina nel vuoto, `1.815` atomi, carica `-5.00 e`).
    - `2myo_wb_ionized.psf` e `2myo_wb_ionized.pdb` (proteina + acqua + ioni, **`21.461` atomi**, `6.535` molecole d'acqua, `23` $\text{Na}^+$, `18` $\text{Cl}^-$, carica totale **`0.0000 e`**).
    - `pbc_box_params.txt` (contiene i vettori della scatola `71.840 × 58.972 × 54.840 Å` e centro `7.735, -8.087, 86.003`).

---

### 4. `scripts/viz_02_solvated_box.tcl` — 📸 Visualizzazione 3D della Scatola Periodica d'Acqua (Slide 3)
- **Come si lancia:**
  ```bash
  vmd -e scripts/viz_02_solvated_box.tcl
  ```
- **Cosa aspettarsi come output VISIVO (Finestra 3D di VMD):**
  - Sfondo bianco con un **parallelepipedo delineato da bordi neri** (`pbc box`).
  - Al centro c'è la **Myotrophin** circondata da un reticolo azzurro trasparente (le `6.535` molecole d'acqua) e da **23 sfere Blu** (ioni Sodio $\text{Na}^+$) e **18 sfere Arancioni** (ioni Cloro $\text{Cl}^-$) distribuite nel volume.

---

### 5. `scripts/03_analyze_wt_baseline.tcl` — Analisi Completa delle Traiettorie MD (`analysisMD.tcl` del Prof. + $Q$-value)
- **Come si lancia:**
  ```bash
  vmd -dispdev text -e scripts/03_analyze_wt_baseline.tcl
  ```
- **Cosa fa:** Integra **tutte le analisi di `/home/olly/UNI/BioMolEng2026/laboratorio/analysisMD.tcl` del Prof. Gautieri** con il calcolo del **$Q$-value di Bekker et al. (2019)**:
  - Allinea tutti i fotogrammi della simulazione sul backbone del fotogramma 0 (`measure fit`).
  - Calcola per il WT a $300\text{ K}$ e a $400\text{ K}$: **RMSD del backbone**, **RMSF dei $C_\alpha$**, **Numero e % di occupazione degli Hydrogen Bonds**, **SASA** (superficie esposta al solvente), **Raggio di Girazione ($R_g$)**, **Percentuale di $\alpha$-elica (`AlphaPerc` tramite STRIDE)**, **$Q$-value nel tempo** (`all-all`, `hydrophilic-all`, `hydrophobic-small`, e separato per **`ANK1`, `ANK2`, `ANK3`, `ANK4`**) e infine il **$\langle Q_i \rangle$ per ciascuno dei 118 amminoacidi** negli ultimi 30% della run a $400\text{ K}$.
- **Cosa aspettarsi come output (File generati in `02_wt_baseline/`):**
  - `rmsd_300K.dat` e `rmsd_400K.dat`: 2 colonne (`Frame`, `RMSD in Å`). A $300\text{ K}$ ci aspettiamo valori bassi (~1.5–2.5 Å); a $400\text{ K}$ valori più alti (~3.0–5.0 Å).
  - `rmsf_300K.dat` e `rmsf_400K.dat`: 2 colonne (`Resid 1..118`, `RMSF in Å`). La curva a $400\text{ K}$ sarà più alta di quella a $300\text{ K}$, con picchi sui loop e nelle zone instabili.
  - `hbonds_300K.dat`, `hbonds_400K.dat` e `hbonds_details_400K.txt`: numero di legami a idrogeno nel tempo e tabella con la `% di esistenza` di ogni singola coppia donatore–accettore!
  - `sasa_300K.dat`, `sasa_400K.dat`, `rgyr_300K.dat`, `rgyr_400K.dat`, `alpha_300K.dat`, `alpha_400K.dat`: evoluzione nel tempo della superficie esposta ($\text{Å}^2$), della compattezza ($R_g \approx 14.5\text{--}15.5\text{ Å}$) e della % di $\alpha$-elica ($\approx 45\text{--}55\%$).
  - `q_time_300K.dat` e `q_time_400K.dat`: 8 colonne (`Frame`, `Q_all`, `Q_phil_all`, `Q_phob_small`, `Q_ANK1`, `Q_ANK2`, `Q_ANK3`, `Q_ANK4`). Parte da $1.000$ al frame 0; a $300\text{ K}$ resta alto ($\sim 0.85\text{--}0.95$), mentre a $400\text{ K}$ scende ($\sim 0.65\text{--}0.80$) mostrando quale dei 4 Ankyrin Repeats cede di più!
  - **`unstable_residues_report.txt`**: La **tabella leggibile a occhio umano** con tutti i 118 amminoacidi, quanti contatti nativi avevano all'inizio, il loro $\langle Q_i \rangle$ finale a $400\text{ K}$ e l'etichetta `>>> UNSTABLE (Q<0.6) <<<`, `WEAK` o `MASKED`!
  - `wt_400K_q_beta.pdb`: L'ultimo fotogramma della proteina a $400\text{ K}$ in cui nella colonna `Beta` di ogni atomo è scritto il valore $\langle Q_i \rangle$ del suo amminoacido (serve per colorare in 3D la proteina nello script 8!).

---

### 6. `scripts/plot_wt_baseline.py` — 📈 Generazione Automatica dei Grafici 2D Vettoriali per le Slide (Slide 4)
- **Come si lancia:**
  ```bash
  python3 scripts/plot_wt_baseline.py
  ```
- **Cosa aspettarsi come output VISIVO (2 File PDF Vettoriali in `06_presentation/`):**
  1. **`06_presentation/fig5_wt_baseline_plots.pdf`**: Pannello vettoriale a 4 quadranti pronto da inserire nella **Slide 4**:
     - *In alto a sinistra:* Curva **RMSD** ($300\text{ K}$ blu vs $400\text{ K}$ rossa).
     - *In alto a destra:* Curva **$Q(t)$** (`all-all` e `hydrophilic-all`).
     - *In basso a sinistra:* Profilo **RMSF** lungo i 118 amminoacidi con le 4 bande verticali colorate di `ANK1`, `ANK2`, `ANK3`, `ANK4`.
     - *In basso a destra:* Profilo **$\langle Q_i \rangle$ a $400\text{ K}$** con la linea tratteggiata rossa a $Q_i = 0.60$, i pallini blu della *Functional Mask* e i **pallini rossi con i nomi dei 5 amminoacidi più instabili dello scaffold**!
  2. **`06_presentation/fig5b_wt_analysisMD_observables.pdf`**: Pannello vettoriale a 4 quadranti con tutte le osservabili di `analysisMD.tcl` del Prof. Gautieri:
     - *In alto a sinistra:* **% di $\alpha$-elica** nel tempo ($300\text{ K}$ vs $400\text{ K}$).
     - *In alto a destra:* **Numero di Hydrogen Bonds** nel tempo ($300\text{ K}$ vs $400\text{ K}$).
     - *In basso a sinistra:* **Raggio di Girazione ($R_g$)** nel tempo.
     - *In basso a destra:* **$Q(t)$ separato per i 4 singoli Ankyrin Repeats (`ANK1` Blu, `ANK2` Ciano, `ANK3` Arancio, `ANK4` Verde)** a $400\text{ K}$!

---

### 7. `scripts/viz_03_wt_instability.tcl` — 📸 Visualizzazione 3D dei Punti Deboli a 400 K (Slide 4)
- **Come si lancia (dopo aver eseguito `03_analyze_wt_baseline.tcl`):**
  ```bash
  vmd -e scripts/viz_03_wt_instability.tcl
  ```
- **Cosa aspettarsi come output VISIVO (Finestra 3D di VMD):**
  - Vedrete due strutture sovrapposte:
    1. **Ombra Grigia Trasparente:** la forma iniziale perfetta della Myotrophin all'istante 0.
    2. **Struttura Opaca Colorata (`Scala Rosso-Bianco-Blu`):** la forma della Myotrophin dopo lo stress termico a $400\text{ K}$, dove il colore di ogni pezzo indica il suo $\langle Q_i \rangle$:
       - **Blu scuro** = Pezzi rimasti super-stabili ($\langle Q_i \rangle \approx 0.90\text{--}1.00$).
       - **Bianco/Rosa** = Pezzi mediamente flessibili ($\langle Q_i \rangle \approx 0.65\text{--}0.75$).
       - **Rosso acceso (con le catene laterali disegnate a bastoncino)** = I **"bulloni instabili" ($\langle Q_i \rangle < 0.65$)** che hanno ceduto al calore e che andremo a rinforzare con le mutazioni nella Fase 3!
