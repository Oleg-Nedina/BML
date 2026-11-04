# 🧬 Computational Design of Thermostable Protein Variants: Myotrophin (`2MYO`)

**Corso:** Biomolecular Modelling Laboratory (2026) — Prof. Alfonso Gautieri (Politecnico di Milano)  
**Team #10:** Oleg Nedina & Alessia Rigoni  
**Proteina Assegnata:** **Myotrophin** (*Rattus norvegicus*) — **PDB ID:** [`2MYO`](https://www.rcsb.org/structure/2MYO) (118 amminoacidi)  
**Presentazione Finale:** 📅 **Mercoledì 16 Dicembre** (durata massima: **10 minuti**)

---

## 🎯 1. Obiettivo del Progetto e Contesto Scientifico

L'obiettivo del progetto è **progettare computazionalmente varianti termostabili** della proteina assegnata (**Myotrophin**, PDB `2MYO`), seguendo il protocollo di riferimento di **Bekker, Ma & Kamiya (2019)** (*"Thermal stability of single-domain antibodies estimated by molecular dynamics simulations"*, Protein Science 28:429–438) ed estendendolo con una **pipeline multi-tool di livello da ricerca (*Research-Grade Workflow*)** ispirata alla `Lecture05` del corso e alla letteratura recente di protein engineering (*FireProt + Rosetta + ProteinMPNN + Consensus Design*).

### L'idea chiave del metodo (Bekker et al., 2019 + Upgrade "Research-Grade")

1. **Perché non misurare direttamente la $T_m$ con la Dinamica Molecolare?**  
   La temperatura di melting ($T_m$) è la temperatura a cui metà della popolazione proteica è folded e metà è unfolded all'equilibrio. Simulare il completo unfolding/folding all'equilibrio vicino alla $T_m$ sperimentale richiederebbe tempi di simulazione inaccessibili (millisecondi).
2. **Simulazioni MD ad alta temperatura + Frazione di Contatti Nativi ($Q$-value):**  
   Eseguendo simulazioni di Dinamica Molecolare (MD) a **temperatura elevata** (es. **400 K**, confrontata con un controllo a **300 K**), si accelera la dinamica conformazionale locale senza distruggere istantaneamente il fold globale (cosa che invece accade a 500 K). Misurando la **frazione di contatti atomici nativi ($Q$-value)** conservati nella parte finale della traiettoria e le fluttuazioni per residuo (**RMSF** e **$Q_i$ per residuo**), è possibile:
   - Identificare i **residui localmente instabili** (i "punti deboli" strutturali che perdono i contatti nativi ad alta temperatura, $\langle Q_i \rangle < 0.6$).
   - Incrociare questi punti deboli dinamici con **strumenti ortogonali di progettazione proteica** (**FoldX**, **Rosetta**, **FireProt 2.0**, **Consenso degli Ankyrin Repeats** e **ProteinMPNN/ThermoMPNN**) per selezionare mutazioni singole e combinate ad alta confidenza.
   - Validare *in silico* le varianti selezionate tramite **simulazioni MD comparative ("matched")** con più repliche indipendenti (*velocity seeds*), dimostrando un aumento statisticamente significativo del **$Q$-value** (globale, per classe chimica e per singolo Ankyrin Repeat) rispetto al Wild-Type (WT).

---

## 🚀 2. Come Trasformare il Progetto in un "Lavoro da Ricerca" (*Research-Grade Upgrades*)

Pur rispettando al 100% la traccia e i vincoli del Prof. Gautieri (FoldX + NAMD + $Q$-value secondo Bekker et al.), il nostro workflow introduce **4 miglioramenti scientifici** che elevano il progetto allo standard di una pubblicazione peer-reviewed:

### 🔬 Upgrade 1: "Orthogonal Consensus Funnel" per la Scelta delle Mutazioni (Fisica + Evoluzione + AI)
Affidarsi a un solo predittore energetico a backbone semi-rigido (come il solo FoldX) espone al rischio di **falsi positivi**. Come mostrato nella `Lecture05` e nel paper del Journal Club presente nel materiale del corso (*ScALDO Engineering, 2025*, che combina **FireProt + Rosetta** per ottenere $\Delta T_m = +6.5^\circ\text{C}$), noi selezioniamo i mutanti facendo convergere **4 pilastri indipendenti**:

| Pilastro | Strumento / Metodo | Cosa valuta (Razionale Biofisico) | Criterio di Selezione |
| :--- | :--- | :--- | :--- |
| **1. Vulnerabilità Dinamica** *(MD-driven)* | **NAMD 3 + VMD (`VMDextensions`)** sul WT a 400 K | Individua i residui che perdono contatti nativi sotto stress termico (approccio *Bekker et al., 2019*) | Residui con $\langle Q_i \rangle_{400\text{K}} < 0.6$ o picchi anomali di RMSF / perdita di H-bond inter-repeat |
| **2. Evoluzione & Architettura a Repeat** *(Phylogeny & Consensus)* | **FireProt 2.0** + **Allineamento Consenso Ankyrin Repeat (`TPLHLAA...`)** | Sfrutta la conservazione evolutiva (*Back-to-Consensus* / *Ancestral Sequence Reconstruction*) e la regolarità strutturale dei 4 Ankyrin Repeats di `2MYO` | Posizioni in cui `2MYO` presenta un residuo "frustrato" che devia dal consenso canonico delle Ankyrin Repeat (es. *DARPins*) |
| **3. Energia Fisica / Empirica** *(Force-Field $\Delta\Delta G$)* | **FoldX 5** (`PositionScan` / `BuildModel`) + **Rosetta** (`Cartesian_ddG` via FireProt/ROSIE/RaSP) | Stima la variazione di energia libera di folding ($\Delta\Delta G_{\text{fold}}$) con funzioni di scoring complementari (empirica in FoldX vs rotamer packing Monte Carlo in Rosetta) | **Consenso energetico:** $\Delta\Delta G_{\text{FoldX}} < -1.0\text{ kcal/mol}$ **E** $\Delta\Delta G_{\text{Rosetta}} < -1.0\text{ kcal/mol}$ |
| **4. Deep Learning Strutturale** *(Inverse Folding & Graph AI)* | **ProteinMPNN / ThermoMPNN** (`Lecture05`) + **DynaMut2** | Valuta la probabilità condizionata della sequenza dato il backbone 3D ($\log P(\text{aa} \mid \text{backbone})$) e l'effetto delle vibrazioni NMA sulla stabilità | $\Delta\log P_{\text{MPNN}} > 0$ e $\Delta\Delta G_{\text{ThermoMPNN/DynaMut2}}$ stabilizzante, con **maschera sui residui funzionali** |

### 🧬 Upgrade 2: Analisi Biofisica "Repeat-Resolved" Specifica per la Myotrophin (`ANK1–ANK4`)
In letteratura biofisica (studi classici del gruppo di **Laura Itzhaki** a Cambridge su Myotrophin `2MYO`, es. *Lowe & Itzhaki, PNAS 2007*), la Myotrophin è il modello archetipico per il folding delle proteine ad Ankyrin Repeat:
- I **repeat C-terminali (Repeat 3 e Repeat 4)** costituiscono il **nucleo di folding primario** (più rigido), mentre i **repeat N-terminali (Repeat 1 e Repeat 2)** sono più plastici.
- La stabilità globale dipende criticamente dalle **interfacce inter-repeat (`ANK1–2`, `ANK2–3`, `ANK3–4`)** che accoppiano cooperativamente i moduli adiacenti.
- **Il nostro valore aggiunto:** Oltre a calcolare il $Q$-value globale e per classe chimica (`all-all`, `hydrophilic-all`, `hydrophobic-small` come in Bekker et al.), calcoleremo il **$Q$-value per singolo Ankyrin Repeat (`Q_ANK1` ... `Q_ANK4`) e per le interfacce inter-repeat**, collegando direttamente le nostre simulazioni MD alle evidenze sperimentali di folding della Myotrophin!

### ⚙️ Upgrade 3: Force Field CHARMM36m e Ciclo Epistatico del Doppio Mutante
1. **Scelta del Force Field (`CHARMM36m`):** Nei file dei tutorial (`/home/olly/UNI/BioMolEng2026/Tutorials/charmm36/`) sono presenti sia il vecchio `CHARMM27` sia il moderno **`CHARMM36m`** (`top_all36_prot.rtf` + `par_all36m_prot.prm`). Useremo **CHARMM36m** perché corregge il noto bias di sovra-stabilizzazione delle $\alpha$-eliche di CHARMM22/27, offrendo una descrizione molto più realistica dell'unfolding parziale di una proteina a solenoide $\alpha$-elica come `2MYO`.
2. **Double-Mutant Cycle (Analisi di Additività vs Epistasi):** Confrontando i due mutanti singoli (`Mut1`, `Mut2`) con il **Doppio Mutante (`Mut1 + Mut2`)**, verificheremo se le due mutazioni agiscono in modo **additivo** ($\Delta Q_{1+2} \approx \Delta Q_1 + \Delta Q_2$) oppure **sinergico/cooperativo (epistasi positiva)** lungo il solenoide degli Ankyrin Repeats, accompagnando i dati con **matrici di differenza delle mappe di contatto ($\Delta C_{ij} = \langle C_{ij}^{\text{Mut}} \rangle - \langle C_{ij}^{\text{WT}} \rangle$)** e test statistici ($t$-test tra le repliche).

---

## 🔍 3. Identikit della Proteina Target: Myotrophin (`2MYO`)

Ecco tutti i dettagli strutturali, biofisici e funzionali su `2MYO` (fondamentali per il **Learning Objective 1: Protein Analysis**):

| Proprietà | Dettaglio per `2MYO` |
| :--- | :--- |
| **Nome Proteina** | **Myotrophin** (nota anche in letteratura come proteina **V-1**) |
| **Organismo** | *Rattus norvegicus* (Ratto) — sequenza identica/altamente conservata nei mammiferi |
| **PDB ID** | [`2MYO`](https://www.rcsb.org/structure/2MYO) (struttura media minimizzata dell'ensemble NMR `1MYO`) |
| **Lunghezza e Catene** | **118 amminoacidi** (~12.5 kDa), **singola catena (`Chain A`)**, stato biologico: **Monomero** |
| **Metodo Sperimentale** | **Solution NMR** (1 modello minimizzato rappresentativo di 90 conformeri; RMSD del backbone ~0.43 Å) |
| **Architettura 3D (Fold)** | **4 Ankyrin-like Repeats (`ANK1`–`ANK4`)** impilati linearmente a formare un solenoide $\alpha$-elica leggermente curvo (~46 Å × 21 Å) |
| **Suddivisione dei 4 Repeat** | **ANK1:** res. 4–32 \| **ANK2:** res. 33–65 \| **ANK3:** res. 66–98 \| **ANK4 (C-cap):** res. 99–118 |
| **Stabilità Sperimentale WT** | $\Delta G_{\text{D-N}}^{\text{H}_2\text{O}} \approx 5.1 - 5.8\text{ kcal/mol}$, $T_m \approx 53 - 56^\circ\text{C}$ (folding cooperativo a 2 stati, *Lowe & Itzhaki*) |

### ⚠️ Due Attenzioni Tecniche Cruciali su `2MYO`

1. **È una struttura NMR (contiene già atomi di idrogeno NMR):**
   - Il file `2MYO.pdb` contiene un singolo modello minimizzato con gli idrogeni NMR e nessuna molecola d'acqua cristallografica. Quando generiamo la topologia CHARMM36m con `psfgen`, dobbiamo selezionare **`protein and noh`** per rimuovere gli idrogeni con nomenclatura NMR e far ricostruire a `psfgen` gli idrogeni con la nomenclatura esatta di CHARMM36m.
   - **Stato delle Istidine e Cisteine in `2MYO`:**
     - **Istidine (`HIS`):** Sono presenti **4 istidine** (**His9, His39, His43, His76**). Negli Ankyrin Repeats l'istidina del motivo conservato `TPLH` (come **His43** in ANK2 e **His76** in ANK3) forma un legame a idrogeno strutturale critico tra il suo anello imidazolico e il backbone del $\beta$-hairpin loop adiacente! In Fase 1 controlleremo in VMD quale azoto ($N_\delta$ o $N_\epsilon$) dona/accetta l'H-bond per assegnare correttamente `HSD` o `HSE` a ciascuna istidina.
     - **Cisteine (`CYS`):** Sono presenti **2 cisteine ridotte** (**Cys10, Cys90**) distanti tra loro (~25 Å): **non formano ponti disolfuro**.
2. **Funzione Biologica e "Functional Masking" (Cosa NON mutare!):**
   - **Ruolo biologico:** La Myotrophin (V-1) regola la polimerizzazione dell'actina ed è implicata nell'ipertrofia cardiaca. Agisce legando e inibendo stericamente l'**Actin-Capping Protein (CP)** (complesso cristallografico nel PDB [`3AAA`](https://www.rcsb.org/structure/3AAA)).
   - **Vincolo di Progettazione (*Functional Masking*):** Dal complesso `3AAA` è noto che V-1 interagisce con la faccia concava di CP tramite i **loop $\beta$-hairpin esposti** e i residui chiave **Trp8, Asp44, Glu78** (oltre ai residui adiacenti dei loop 5–12, 42–47, 75–81). **Tutti questi residui funzionali saranno "mascherati" (esclusi dalle mutazioni)**: interverremo esclusivamente sullo scaffold strutturale delle eliche antiparallele e dei contatti inter-repeat per aumentare la termostabilità preservando intatta l'attività biologica!

---

## 🛒 4. "Lista della Spesa" (Checklist Operativa del Progetto)

### Punto 0 — Setup del Sistema e dell'Ambiente di Lavoro (`00_setup`)
- [ ] **0.1** Verificare il funzionamento di **VMD 2.0.0** (`/usr/local/bin/vmd`), **NAMD 3.0.3** (`/usr/local/bin/namd3`) e il caricamento automatico di **`VMDextensions.tcl`** da `~/.vmdrc`.
- [ ] **0.2** Recuperare e configurare l'eseguibile **`foldx`** (estraendolo dalla VM `.vdi` del corso o scaricando la versione Academic da `foldxsuite.crg.eu`) e renderlo eseguibile nel `PATH`.
- [ ] **0.3** Configurare l'ambiente Python locale (`numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`) per le analisi statistiche e i grafici *publication-grade*.
- [ ] **0.4** Verificare l'accesso ai tool complementari di consensus design (**FireProt 2.0**, **Rosetta / RaSP / PremPS**, **ProteinMPNN / ThermoMPNN** via Colab/Web Server, **DynaMut2**).

### Fase 1 — Studio della Proteina, Architettura Ankyrin e Setup Topologico (`01_structures/`)
- [ ] **1.1** Scaricare [`2MYO.pdb`](https://files.rcsb.org/download/2MYO.pdb) (e [`3AAA.pdb`](https://files.rcsb.org/download/3AAA.pdb) come riferimento per l'interfaccia funzionale) in `01_structures/raw/`.
- [ ] **1.2** Ispezionare `2MYO` in VMD: mappare i 4 Ankyrin Repeats (`ANK1`–`ANK4`), verificare lo stato di protonazione specifico (`HSD` vs `HSE`) di **His9, His39, His43, His76** e definire la lista dei residui funzionali da proteggere (*Functional Mask*).
- [ ] **1.3** Eseguire l'**allineamento di sequenza dei 4 Ankyrin Repeats di `2MYO`** contro il motivo consenso universale degli Ankyrin Repeats (`G-X-TPLHLAA...`) per individuare le posizioni "frustrate" non canoniche.
- [ ] **1.4** Generare la topologia **CHARMM36m** del Wild-Type (`2myo_wt.psf` e `2myo_wt.pdb`) con `psfgen`, solvatarla in una scatola d'acqua esplicita TIP3P (padding 12 Å) e ionizzarla a **0.15 M NaCl** (`autoionize`).

### Fase 2 — Benchmark Computazionale e Simulazione Baseline del WT (`02_wt_baseline/`)
- [ ] **2.1** Eseguire il **benchmark di velocità** con **NAMD 3** (`namd3 +p<N_cores>`) sul sistema solvatato e misurare i **`ns/day`**.
- [ ] **2.2** Definire e documentare il **Computational Budget** ($\text{Sistemi} \times \text{Temperature} \times \text{Repliche} \times \text{Durata}$).
- [ ] **2.3** Eseguire le simulazioni baseline del **WT a 300 K** (controllo fisiologico) e ad **alta temperatura (400 K, più eventuale check a 360–380 K)** per indurre stress termico locale senza denaturazione istantanea.
- [ ] **2.4** Analizzare il WT ad alta temperatura:
  - Ricostruire PBC (`pbc wrap`) e allineare i frame sul frame iniziale.
  - Calcolare **RMSD**, **RMSF per residuo**, **$Q$-value per singolo Repeat (`ANK1`–`ANK4`)** e **$Q$-value per singolo residuo ($Q_i$)** per identificare i residui localmente instabili ($\langle Q_i \rangle < 0.6$).

### Fase 3 — Multi-Tool Consensus Screening & Modellazione dei Mutanti (`03_foldx_screening/`)
- [ ] **3.1** **Screening Evolutivo & Ibrido (FireProt 2.0 + Ankyrin Consensus):** Sottomettere `2MYO` a **FireProt 2.0** (ramo evolutivo + energetico FoldX/Rosetta) e incrociare l'output con l'allineamento consenso dei 4 Ankyrin Repeats.
- [ ] **3.2** **Screening AI Strutturale (ProteinMPNN / ThermoMPNN & DynaMut2):** Valutare le posizioni candidate (escludendo i residui funzionali mascherati Trp8, Asp44, Glu78) con ThermoMPNN/ProteinMPNN e DynaMut2.
- [ ] **3.3** **Screening Sistematico FoldX (Tool Ufficiale del Corso):**
  - Eseguire `RepairPDB` su `2MYO.pdb` $\rightarrow$ `2MYO_Repair.pdb`.
  - Eseguire `PositionScan` sulle posizioni instabili/candidate e raccogliere i profili completi di $\Delta\Delta G_{\text{FoldX}}$.
- [ ] **3.4** **Tabella di Consenso Multi-Tool & Costruzione 3D (`BuildModel`):**
  - Costruire la matrice decisionale che incrocia **Instabilità MD ($Q_i$)** + **FoldX $\Delta\Delta G$** + **Rosetta/FireProt $\Delta\Delta G$** + **ThermoMPNN/DynaMut2** + **Ankyrin Consensus**.
  - Selezionare le **2 migliori mutazioni singole (`Mut1`, `Mut2`)** e il **Doppio Mutante combinato (`Mut1 + Mut2`)**.
  - Generare le strutture 3D con `BuildModel` di **FoldX** (`numberOfRuns=3..5`) ed estrarre la scomposizione energetica (`Dif_2MYO_Repair.fxout`: H-bond, VdW, Electrostatics, Solvation, Entropy).

### Fase 4 — Simulazioni MD Comparative ("Matched Comparison", `04_md_simulations/`)
- [ ] **4.1** Preparare **WT**, **Mut1**, **Mut2** e **Doppio Mutante (`Mut1+Mut2`)** con **esattamente lo stesso protocollo CHARMM36m + TIP3P + 0.15 M NaCl**.
- [ ] **4.2** Eseguire Minimizzazione + Equilibrazione NPT + Produzione MD ad alta temperatura per tutti i 4 sistemi con **almeno 3 repliche indipendenti (`seed` diversi)** per ciascun sistema.

### Fase 5 — Analisi Statistica Avanzata, Ciclo Epistatico e Presentazione (`05_analysis/` & `06_presentation/`)
- [ ] **5.1** Calcolare il **$Q$-value** nel tempo e mediarlo sulla finestra finale equilibrata per tutte le repliche (**Media $\pm$ SD** e **$p$-value**):
  - $Q$ globale e per classi chimico-fisiche (**`all-all`**, **`hydrophilic-all`**, **`hydrophobic-small`** — *Bekker et al. Table II*).
  - $Q$ risolto per dominio (**`ANK1`**, **`ANK2`**, **`ANK3`**, **`ANK4`** e interfacce inter-repeat).
- [ ] **5.2** Analizzare l'**Additività / Epistasi del Doppio Mutante** ($\Delta Q_{1+2}$ vs $\Delta Q_1 + \Delta Q_2$) e contare il numero di residui instabili residui ($\langle Q_i \rangle < 0.6$) in ciascuna variante.
- [ ] **5.3** Calcolare le osservabili meccaniche/strutturali: **RMSF comparativo**, **Matrici di differenza dei contatti ($\Delta C_{ij}$)**, **Occupancy di H-bond/Salt Bridges** e conservazione della **frazione di $\alpha$-elica (`helicity`)**.
- [ ] **5.4** Realizzare le slide per il **16 Dicembre** (max 10 minuti) con divisione chiara dei contributi di Oleg e Alessia.

---

## 📂 5. Struttura della Repository `BML` e Divisione del Lavoro

```text
BML/
├── README.md                  # Guida completa, Architettura Research-Grade, Checklist e Roadmap
├── .gitignore                 # Esclude file pesanti (.dcd, .coor, .vel, .log, binari) da Git
├── 01_structures/             # Struttura PDB 2MYO (e 3AAA), studio Ankyrin consensus, PSF/PDB CHARMM36m
├── 02_wt_baseline/            # Benchmark NAMD (ns/day), MD WT (300K vs 400K), mappa Q_i e RMSF
├── 03_foldx_screening/        # Output FoldX (RepairPDB, PositionScan, BuildModel) + FireProt/Rosetta/MPNN
├── 04_md_simulations/         # Configurazioni NAMD (.conf) per il confronto WT vs Mut1, Mut2, Double (multi-seed)
├── 05_analysis/               # Risultati Q-value (globale, chimico, per-repeat), ΔC_ij, grafici e test statistici
├── 06_presentation/           # Slide, render 3D VMD ad alta risoluzione e scaletta per il 16 Dicembre
└── scripts/                   # Script TCL e Python riutilizzabili e commentati
```

---

## 🛠️ 6. Mappa Completa dei Tool (Locali + Research Web Servers)

### A. Software Installati Localmente e File del Corso
| Tool / Risorsa | Percorso sul Sistema | Ruolo nel Progetto |
| :--- | :--- | :--- |
| **NAMD 3.0.3** (Multicore) | `/usr/local/bin/namd3` | Motore MD per minimizzazione, equilibrazione e produzione |
| **VMD 2.0.0** | `/usr/local/bin/vmd` | Visualizzazione 3D, `psfgen`, `solvate`, `autoionize` e analisi TCL |
| **VMDextensions** (T. Giorgino) | `/usr/local/lib/vmd/plugins/VMDextensions.tcl` | Funzioni `prepareNativeContacts`, `measureNativeContacts`, `forFrames`, `helicity` |
| **Force Field CHARMM36m** | `/home/olly/UNI/BioMolEng2026/Tutorials/charmm36/` | Topologia `top_all36_prot.rtf` e parametri `par_all36m_prot.prm` |
| **FoldX Suite** | Da installare in `/usr/local/bin/foldx` (vedi Punto 0) | `RepairPDB`, `PositionScan` e `BuildModel` per $\Delta\Delta G$ e costruzione 3D mutanti |
| **Tutorial Native Contacts ($Q$)** | `/home/olly/UNI/BioMolEng2026/Tutorials/nativecontacts/` | Script di riferimento `native_contacts.tcl` |
| **Tutorial FoldX** | `/home/olly/UNI/BioMolEng2026/Tutorials/foldx_test/` | Esempi di configurazione per `PositionScan`, `BuildModel` e `Stability` |

### B. Strumenti "Research-Grade" Complementari (Gratuiti via Web Server / Colab)
| Tool | URL / Accesso | Perché usarlo insieme a FoldX e MD |
| :--- | :--- | :--- |
| **FireProt 2.0** | [loschmidt.chemi.muni.cz/fireprot](https://loschmidt.chemi.muni.cz/fireprot/) | **Citato nel paper ScALDO 2025 del corso!** Combina automaticamente *Back-to-Consensus* evolutivo con filtri energetici **FoldX** e **Rosetta** |
| **Rosetta ($\Delta\Delta G$ / RaSP)** | Integrato in FireProt / [ROSIE](https://rosie.rosettacommons.org/) / Colab **RaSP** (*Rapid Stability Prediction*) | Valuta il *rotamer packing* con l'energia fisica/statistica di **Rosetta** (`Lecture05`, Slide 48–59), ortogonale a FoldX |
| **ThermoMPNN / ProteinMPNN** | Google Colab ufficiale ([ThermoMPNN](https://colab.research.google.com/github/Kuhlman-Lab/ThermoMPNN/blob/main/ThermoMPNN.ipynb)) | Rete neurale Graph-based vista nella `Lecture05` (Slide 65–69), addestrata su megascala sperimentale di $\Delta\Delta G$ di termostabilità |
| **DynaMut2 / DDMut** | [biosig.lab.uq.edu.au/dynamut2](https://biosig.lab.uq.edu.au/dynamut2/) | Combina grafi strutturali e **Normal Mode Analysis (NMA)** per stimare l'effetto della mutazione sia su $\Delta\Delta G$ sia sulla flessibilità |

---

## 🗺️ 7. Roadmap Operativa Passo-Passo (Dal Punto 0 alla Presentazione)

---

### 📍 PUNTO 0: Setup del Sistema e dell'Ambiente di Lavoro

Prima di iniziare la Fase 1, verifichiamo che tutta la toolchain locale sia pronta e funzionante.

#### 0.1 Verifica di VMD, NAMD 3 e `VMDextensions.tcl`
Sul vostro sistema Linux (come documentato in `/home/olly/UNI/BioMolEng2026/dependencies.md`):
- **NAMD 3.0.3** è già installato in `/usr/local/bin/namd3`.
- **VMD 2.0.0** è già installato in `/usr/local/bin/vmd`.
- La libreria **`VMDextensions.tcl`** è presente in `/usr/local/lib/vmd/plugins/VMDextensions.tcl` e viene caricata automaticamente all'avvio di VMD tramite il file `~/.vmdrc`.
- **Test rapido da terminale:**
  ```bash
  namd3 | head -n 5
  vmd -dispdev text -eofexit <<< "info procs prepareNativeContacts"
  ```
  Se il secondo comando risponde `prepareNativeContacts`, VMD e le estensioni per il calcolo del $Q$-value sono pronte!

#### 0.2 Installazione del Binario `foldx`
Nei file di configurazione dei tutorial (`/home/olly/UNI/BioMolEng2026/Tutorials/foldx_test/build_model/config.cfg`), il binario di FoldX nella macchina virtuale (`.vdi`) del corso si trovava in `/home/biomolecular/Software/foldx`. Per averlo nativamente sul vostro Linux avete due opzioni semplicissime:
- **Opzione 1 (Estrazione dalla VM `.vdi` del corso, come fatto in `dependencies.md`):**
  Se avete ancora il file `.vdi` della VM del corso, potete montarlo con `guestmount` e copiare `/mnt/vdi/home/biomolecular/Software/foldx` (e `rotabase.txt`, se usate FoldX 4) in `/usr/local/bin/`:
  ```bash
  sudo cp /mnt/vdi/home/biomolecular/Software/foldx /usr/local/bin/foldx
  sudo chmod +x /usr/local/bin/foldx
  ```
- **Opzione 2 (Download diretto da FoldX Suite — Licenza Accademica Gratuita):**
  Registrarsi con la mail istituzionale `@mail.polimi.it` su [https://foldxsuite.crg.eu/](https://foldxsuite.crg.eu/), scaricare il pacchetto Linux (`foldx5Linux64.tar.gz`) ed estrarre il binario `foldx` in `/usr/local/bin/foldx`.

#### 0.3 Ambiente Python per Analisi Statistica e Grafici da Pubblicazione
Per elaborare gli output `.dat` di VMD/NAMD e produrre figure di qualità da rivista scientifica (curve $Q(t)$ con bande di deviazione standard $\mu \pm \sigma$, heatmap $\Delta C_{ij}$, barplot con $p$-value), assicuratevi di avere disponibili:
```bash
python3 -c "import numpy, scipy, pandas, matplotlib; print('Python OK!')"
```

---

### 📍 FASE 1: Analisi della Proteina Target, Ankyrin Consensus e Setup Topologico (`01_structures/`)

#### 1.1 Download dei PDB (`2MYO` e il complesso funzionale `3AAA`)
Scarichiamo sia la struttura NMR della Myotrophin (`2MYO.pdb`) sia la struttura cristallografica del complesso Myotrophin–Actin Capping Protein (`3AAA.pdb`, utile per visualizzare e mascherare l'interfaccia funzionale):
```bash
wget https://files.rcsb.org/download/2MYO.pdb -O 01_structures/2MYO.pdb
wget https://files.rcsb.org/download/3AAA.pdb -O 01_structures/3AAA_complex_ref.pdb
```

#### 1.2 Analisi Strutturale dei 4 Ankyrin Repeats e Stato delle 4 Istidine
Aprendo `01_structures/2MYO.pdb` in VMD:
1. **Architettura dei 4 Ankyrin Repeats (`ANK1`–`ANK4`):**
   Ogni Ankyrin Repeat (~33 residui) è formato da **due $\alpha$-eliche antiparallele** (elica interna concava ed elica esterna convessa) collegate da un loop corto e seguite da un lungo **loop $\beta$-hairpin/turn** che proietta verso l'esterno (dove lega la Capping Protein).
   - **ANK1 (N-cap):** residui `4–32` (Elica 1A: `7–14`, Elica 1B: `17–24`)
   - **ANK2 (Core):** residui `33–65` (Motivo `TPLH`: `36–39` / `40–43`, Elica 2A: `39–46`, Elica 2B: `50–58`)
   - **ANK3 (Core):** residui `66–98` (Motivo `TPLH`: `73–76`, Elica 3A: `72–79`, Elica 3B: `83–91`)
   - **ANK4 (C-cap):** residui `99–118` (Elica 4A: `105–112`, C-tail: `113–118`)
2. **Controllo della Protonazione delle 4 Istidine (`His9`, `His39`, `His43`, `His76`):**
   Invece di usare un alias cieco uguale per tutte le istidine, controlliamo in VMD i legami a idrogeno di ciascuna istidina (`HSD` = protone su `ND1`; `HSE` = protone su `NE2`) per preservare la rete nativa di H-bond degli Ankyrin Repeats!
3. **Mappa del "Functional Masking" (Residui da NON mutare):**
   Allineando in VMD `2MYO.pdb` sulla catena V-1 di `3AAA_complex_ref.pdb`, identifichiamo tutti i residui a distanza $< 4.5\text{ Å}$ dall'Actin-Capping Protein (in particolare i loop $\beta$-hairpin e la triade **Trp8, Asp44, Glu78**). Questi residui vengono inseriti nella nostra *Exclusion Mask*.

#### 1.3 Allineamento dei 4 Repeat di `2MYO` contro il Consenso Universale Ankyrin
In letteratura (*Mosavi et al., PNAS 2002*; *Binz & Plückthun, DARPins*), il motivo consenso di un Ankyrin Repeat ottimizzato per la massima stabilità termica è:
$$\text{Consenso ANK (33 aa):} \quad \texttt{G-X-T-P-L-H-L-A-A-X-X-G-H-L-E-I-V-E-V-L-L-K-X-G-A-D-V-N-A-X-D-X-X}$$
Allineando le sequenze di `ANK1`, `ANK2`, `ANK3` e `ANK4` di `2MYO` contro questo consenso, individuiamo subito dove la Myotrophin naturale presenta residui "sub-ottimali" (es. cariche non compensate, assenza di Prolina/Treonina nel motivo `TPLH`, o impaccamento idrofobico incompleto tra le eliche antiparallele).

#### 1.4 Generazione Topologia `CHARMM36m` e Solvatazione Esplicita (`psfgen` + `solvate` + `autoionize`)
Usiamo il force field **CHARMM36m** (`/home/olly/UNI/BioMolEng2026/Tutorials/charmm36/top_all36_prot.rtf` e `par_all36m_prot.prm`), rimuovendo gli idrogeni NMR (`protein and noh`), solvatando in un box TIP3P con padding **12 Å** e aggiungendo **0.15 M NaCl**.

---

### 📍 FASE 2: Benchmark, Computational Budget e Vulnerabilità Dinamica del WT (`02_wt_baseline/`)

#### 2.1 Benchmark di Velocità (`ns/day`) e Computational Budget
Lanciamo una run di test di 5.000 step (10 ps) con `namd3 +p<N>` sul sistema solvatato per misurare i **`ns/day`** reali del nostro hardware e compilare la tabella del **Computational Budget** (Slide 21 di `Laboratory_Project_2026.pdf`):
$$\text{Total MD Time} = N_{\text{systems}} \times N_{\text{temperatures}} \times N_{\text{replicas}} \times t_{\text{sim}}$$

#### 2.2 Simulazione Baseline del WT (300 K vs 400 K) e Mappatura dei Residui Instabili ($Q_i$)
1. Eseguiamo il WT a **300 K** (riferimento fisiologico stabile) e a **400 K** (stress termico accelerato, come in *Bekker et al., 2019*).
2. Sulla traiettoria a 400 K (dopo `pbc wrap` e allineamento sui $C_\alpha$ del frame 0) calcoliamo:
   - **RMSD globale** e **RMSD/Q-value separato per ciascuno dei 4 Ankyrin Repeats (`ANK1`–`ANK4`)** $\rightarrow$ Verifichiamo se i repeat N-terminali o C-terminali cedono per primi!
   - **RMSF per residuo ($C_\alpha$)**.
   - **$Q$-value per singolo residuo ($\langle Q_i \rangle$)** negli ultimi nanosecondi della traiettoria (con cutoff $4.5\text{ Å}$ sugli atomi pesanti `noh`, seguendo esattamente Bekker et al.): tutti i residui non terminali e non funzionali con **$\langle Q_i \rangle < 0.6$** (o forte calo rispetto a 300 K) diventano i nostri **target primari di stabilizzazione**!

---

### 📍 FASE 3: Multi-Tool Consensus Screening & Costruzione dei Mutanti (`03_foldx_screening/`)

Questa è la fase in cui il nostro approccio "da ricerca" brilla maggiormente:

1. **Step 3.1 — Screening con FoldX (Baseline Obbligatoria del Corso):**
   - `foldx --command=RepairPDB --pdb=2MYO.pdb` $\rightarrow$ genera `2MYO_Repair.pdb`.
   - `foldx --command=PositionScan --pdb=2MYO_Repair.pdb --positions=...` sulle posizioni instabili emerse dalla MD e dall'analisi del consenso Ankyrin (o sull'intero scaffold non funzionale).
2. **Step 3.2 — Screening Ortogonale con FireProt 2.0, Rosetta e ThermoMPNN/DynaMut2:**
   - Sottomettiamo `2MYO` a **FireProt 2.0** (che esegue automaticamente l'analisi filogenetica *Back-to-Consensus* + calcoli energetici **FoldX** e **Rosetta**).
   - Eseguiamo **ThermoMPNN / ProteinMPNN** e **DynaMut2** su `2MYO_Repair.pdb`.
3. **Step 3.3 — Matrice di Convergenza ("Orthogonal Funnel") e Scelta dei Candidati:**
   Selezioniamo **2 mutazioni singole (`Mut1`, `Mut2`)** (idealmente situate in punti chiave dei repeat/inter-repeat che mostrano bassa stabilità dinamica $Q_i$ nella MD del WT) e il corrispondente **Doppio Mutante (`Mut1 + Mut2`)** che soddisfano contemporaneamente:
   - Esclusione dai siti di legame con Capping Protein (*Functional Masking* rispettato).
   - $\Delta\Delta G_{\text{FoldX}} < -1.0\text{ kcal/mol}$.
   - Concordanza di segno ($\Delta\Delta G < 0$) con **Rosetta / FireProt / ThermoMPNN / DynaMut2** e/o riallineamento verso il **Consenso Ankyrin Repeat**.
4. **Step 3.4 — Generazione dei Modelli 3D con FoldX `BuildModel`:**
   - Creiamo `individual_list.txt` con le 2 mutazioni singole e la doppia mutazione (ricordando il `;` finale su ogni riga!).
   - Eseguiamo `BuildModel` (`numberOfRuns=3`) ed estraiamo da `Dif_2MYO_Repair.fxout` la scomposizione energetica dettagliata ($\Delta E_{\text{H-bond}}$, $\Delta E_{\text{vdW}}$, $\Delta E_{\text{electrostatics}}$, $\Delta E_{\text{solvation}}$) per spiegare chimicamente perché ogni mutazione stabilizza il fold.

---

### 📍 FASE 4: Simulazioni MD Comparative ("Matched Comparison", `04_md_simulations/`)

1. **Preparazione Identica ("Matched Protocol"):**
   Costruiamo i sistemi solvatati e ionizzati (CHARMM36m + TIP3P + 0.15 M NaCl) per i **4 sistemi**:
   - **`WT`** (`2MYO` Wild-Type)
   - **`Mut1`** (Mutante Singolo 1)
   - **`Mut2`** (Mutante Singolo 2)
   - **`Double`** (Doppio Mutante `Mut1 + Mut2`)
2. **Repliche Indipendenti (`Multi-Seed`):**
   Per ciascuno dei 4 sistemi eseguiamo **almeno 3 repliche indipendenti** (variando `seed` in NAMD) alla temperatura di stress termico scelta (es. 400 K), con identico protocollo di minimizzazione, equilibrazione NPT a 300 K e produzione a 400 K.

---

### 📍 FASE 5: Analisi Statistica, Meccanismo Fisico e Presentazione (`05_analysis/` & `06_presentation/`)

1. **Tabella Globale e per Classi Chimiche (*Bekker et al. Table II* + Upgrade Multi-Tool):**

| Variante | $\Delta\Delta G_{\text{FoldX}}$ | $\Delta\Delta G_{\text{Rosetta/MPNN}}$ | RMSD finale (Å) | $Q$: `all-all` | $Q$: `hydrophilic-all` | $Q$: `hydrophobic-small` | N° Residui Instabili ($Q_i < 0.6$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **WT (`2MYO`)** | $0.00$ | $0.00$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $N_{\text{WT}}$ |
| **Mutante 1** | $< 0$ | $< 0$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $< N_{\text{WT}}$ |
| **Mutante 2** | $< 0$ | $< 0$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $< N_{\text{WT}}$ |
| **Doppio (`1+2`)** | $< 0$ | $< 0$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\ll N_{\text{WT}}$ |

2. **Analisi "Repeat-Resolved" (`Q_ANK1` ... `Q_ANK4`) e Additività/Epistasi:**
   - Mostriamo quali Ankyrin Repeats vengono stabilizzati direttamente e se la stabilizzazione si propaga cooperativamente ai repeat adiacenti (effetto allosterico inter-repeat!).
   - Calcoliamo il ciclo epistatico del doppio mutante e i $p$-value tra WT e mutanti.
3. **Evidenze Strutturali di Supporto:**
   - ** Profilo $\Delta\text{RMSF}$** ($\text{RMSF}_{\text{Mut}} - \text{RMSF}_{\text{WT}}$) lungo i 118 residui.
   - **Matrice di Differenza dei Contatti ($\Delta C_{ij}$)** tra mutante e WT.
   - **Occupancy (%) dei nuovi Hydrogen Bonds / Salt Bridges** e conservazione della **frazione di $\alpha$-elica (`helicity`)**.

---

## 📚 8. Riferimenti Bibliografici Chiave (Da Citare nelle Slide!)

1. **Protocollo MD & $Q$-value (Paper di Riferimento del Corso):**  
   Bekker GJ, Ma B, Kamiya N (2019). *Thermal stability of single-domain antibodies estimated by molecular dynamics simulations.* **Protein Science**, 28(2):429–438.
2. **Struttura NMR della Myotrophin (`2MYO`):**  
   Yang Y, Nanduri S, Sen S, Qin J (1998). *The structural basis of ankyrin-like repeat function as revealed by the solution structure of myotrophin.* **Structure**, 6(5):619–626.
3. **Meccanismo di Folding e Stabilità degli Ankyrin Repeats di Myotrophin (`2MYO`):**  
   Lowe AR, Itzhaki LS (2007). *Rational redesign of the folding pathway of a modular protein.* **PNAS**, 104(8):2679–2684.  
   Mosavi LK, Minor DL, Peng ZY (2002). *Consensus-derived structural determinants of the ankyrin repeat motif.* **PNAS**, 99(25):16029–16034.
4. **Complesso Funzionale Myotrophin / Actin-Capping Protein (`3AAA`):**  
   Takeda S et al. (2010). *Two distinct mechanisms for actin capping protein regulation—steric and allosteric inhibition.* **PLoS Biology**, 8(7):e1000416.
5. **Approccio Combinato FireProt + Rosetta + FoldX (Dal Materiale del Corso 2025):**  
   Bednar D et al. (2015) / Musil M et al. (2024). *FireProt 2.0: web server for automated design of thermostable proteins.* **Briefings in Bioinformatics** / **Nucleic Acids Research**.
6. **Deep Learning per Stabilità Proteica (Dalla `Lecture05` del Corso):**  
   Dauparas J et al. (2022). *Robust deep learning–based protein sequence design using ProteinMPNN.* **Science**, 378:49–56.  
   Dieckhaus H, Brocidiacono M, Randolph NZ, Kuhlman B (2024). *Transfer learning to leverage larger datasets for improved prediction of protein stability changes (ThermoMPNN).* **PNAS**, 121(6):e2314853121.
