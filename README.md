# 🧬 Computational Design of Thermostable Protein Variants: Myotrophin (`2MYO`)

**Corso:** Biomolecular Modelling Laboratory (2026) — Prof. Alfonso Gautieri (Politecnico di Milano)  
**Team #10:** Oleg Nedina & Alessia Rigoni  
**Proteina Assegnata:** **Myotrophin** (*Rattus norvegicus*) — **PDB ID:** [`2MYO`](https://www.rcsb.org/structure/2MYO) (118 amminoacidi)  
**Presentazione Finale:** 📅 **Mercoledì 16 Dicembre** (durata massima: **10 minuti**)

---

## 🎯 1. Obiettivo del Progetto e Contesto Scientifico

L'obiettivo del progetto è **progettare computazionalmente varianti termostabili** della proteina assegnata (**Myotrophin**, PDB `2MYO`), seguendo e adattando il protocollo pubblicato nello studio di riferimento di **Bekker, Ma & Kamiya (2019)** (*"Thermal stability of single-domain antibodies estimated by molecular dynamics simulations"*, Protein Science 28:429–438).

### L'idea chiave del metodo (Bekker et al., 2019 + Workflow del Corso)
1. **Perché non misurare direttamente la $T_m$ con la Dinamica Molecolare?**  
   La temperatura di melting ($T_m$) è la temperatura a cui metà della popolazione proteica è folded e metà è unfolded all'equilibrio. Simulare il completo unfolding/folding all'equilibrio vicino alla $T_m$ sperimentale richiederebbe tempi di simulazione inaccessibili (millisecondi).
2. **Simulazioni MD ad alta temperatura + Frazione di Contatti Nativi ($Q$-value):**  
   Eseguendo simulazioni di Dinamica Molecolare (MD) a **temperatura elevata** (es. **400 K**, confrontata con un controllo a **300 K**), si accelera la dinamica conformazionale locale senza distruggere istantaneamente il fold globale (cosa che invece accade a 500 K). Misurando la **frazione di contatti atomici nativi ($Q$-value)** conservati nella parte finale della traiettoria e le fluttuazioni per residuo (**RMSF** e **$Q_i$ per residuo**), è possibile:
   - Identificare i **residui localmente instabili** (i "punti deboli" strutturali che perdono i contatti nativi ad alta temperatura).
   - Formulare **ipotesi di mutazione puntiforme** (singole e combinate/doppie) che migliorino le interazioni locali (es. rimozione di repulsioni elettrostatiche, creazione di ponti salini/H-bond, migliore impaccamento idrofobico).
   - Effettuare uno **screening energetico rapido con FoldX** ($\Delta\Delta G = \Delta G_{\text{mut}} - \Delta G_{\text{WT}} < 0\text{ kcal/mol}$) per selezionare le mutazioni più promettenti.
   - Validare *in silico* le varianti selezionate tramite **simulazioni MD comparative ("matched")** con più repliche indipendenti (*velocity seeds*), dimostrando un aumento del **$Q$-value** (e una riduzione dell'instabilità locale) rispetto al Wild-Type (WT).

---

## 🔍 2. Identikit della Proteina Target: Myotrophin (`2MYO`)

Prima di toccare qualsiasi software, ecco tutto ciò che dovete sapere sulla vostra proteina (fondamentale per l'**Introduzione** della presentazione e per il **Learning Objective 1: Protein Analysis**):

| Proprietà | Dettaglio per `2MYO` |
| :--- | :--- |
| **Nome Proteina** | **Myotrophin** (nota anche in letteratura come proteina **V-1**) |
| **Organismo** | *Rattus norvegicus* (Ratto) |
| **PDB ID** | [`2MYO`](https://www.rcsb.org/structure/2MYO) (collegata all'ensemble NMR `1MYO`) |
| **Lunghezza e Catene** | **118 amminoacidi** (~12.5 kDa), **singola catena (`Chain A`)**, stato biologico: **Monomero** |
| **Metodo Sperimentale** | **Solution NMR** (struttura media minimizzata rappresentativa di 90 conformeri calcolati; precisione del backbone ~0.43 Å) |
| **Architettura 3D (Fold)** | **4 Ankyrin-like Repeats (ANK)** impilati linearmente a formare un solenoide $\alpha$-elica leggermente curvo (~46 Å di lunghezza × 21 Å di larghezza) |
| **Meccanismo di Folding** | Modello biofisico classico per le proteine a repeat: segue un **folding cooperativo a due stati** ($\text{Unfolded} \rightleftharpoons \text{Folded}$) |

### ⚠️ Due Attenzioni Tecniche Cruciali su `2MYO`
1. **È una struttura NMR (non X-ray!):**
   - A differenza delle strutture cristallografiche a raggi X, il file `2MYO.pdb` **contiene già gli atomi di idrogeno** determinati/modellati via NMR e non contiene acqua cristallografica. Quando userete `psfgen` per generare la topologia CHARMM, è buona pratica selezionare solo gli atomi pesanti (`protein and noh`) oppure lasciare che `psfgen` rigeneri gli idrogeni secondo la nomenclatura standard del force field CHARMM.
   - Verificate sempre che nel file scaricato ci sia un singolo modello (in `2MYO` c'è la struttura media minimizzata, quindi 1 modello, a differenza di `1MYO` che contiene l'ensemble).
2. **Funzione Biologica e Vincoli Funzionali (Cosa NON mutare!):**
   - **Ruolo biologico:** La Myotrophin (V-1) è una proteina citosolica implicata nell'innesco dell'**ipertrofia cardiaca** e nella regolazione della polimerizzazione dell'actina. Agisce legando e inibendo stericamente l'**Actin-Capping Protein (CP)** (potete vedere il complesso cristallografico V-1/CP nel PDB [`3AAA`](https://www.rcsb.org/structure/3AAA)), impedendo a CP di bloccare l'estremità *barbed* dei filamenti di actina. Mostra inoltre omologia strutturale con i domini ankyrin di $\text{I}\kappa\text{B}\alpha$ (interazione con proteine NF-$\kappa$B).
   - **Regola d'oro per il Mutation Design:** L'interazione tra Myotrophin e Capping Protein avviene tramite i **loop estesi** che sporgono dallo scaffold delle $\alpha$-eliche degli ankyrin repeats, coinvolgendo residui chiave d'interfaccia (tra cui **Trp8**, **Asp44**, **Glu78**). **Non dovete mutare i residui essenziali per il binding funzionale!** Le mutazioni stabilizzanti vanno progettate sullo scaffold degli ankyrin repeats (eliche antiparallele e contatti inter-repeat lontano dall'interfaccia attiva) per aumentare la stabilità termica senza compromettere la funzione biologica.

---

## 🛒 3. "Lista della Spesa" (Checklist Operativa del Progetto)

Usate questa lista come traccia operativa da spuntare passo dopo passo:

### Fase 1 — Studio della Proteina e Setup Iniziale
- [ ] **1.1** Scaricare [`2MYO.pdb`](https://files.rcsb.org/download/2MYO.pdb) da RCSB PDB e salvarlo in `01_structures/raw/`.
- [ ] **1.2** Studiare la struttura in **VMD**: identificare i 4 Ankyrin Repeats (eliche antiparallele, $\beta$-turn/loop), N-term e C-term, verificare che non ci siano residui o atomi pesanti mancanti.
- [ ] **1.3** Controllare lo stato di protonazione delle **Istidine (`HIS`)** (scegliere `HSD`, `HSE` o `HSP` in base ai legami a idrogeno locali) e verificare eventuali **Cisteine (`CYS`)**.
- [ ] **1.4** Mappare i residui funzionali (interfaccia di legame con Capping Protein, es. Trp8, Asp44, Glu78 e loop esposti) da **escludere** dalle mutazioni.
- [ ] **1.5** Generare la topologia CHARMM del Wild-Type (`2myo_wt.psf` e `2myo_wt.pdb`) tramite `psfgen` in VMD e solvatarla/ionizzarla secondo il protocollo scelto (Box d'acqua esplicita TIP3P + 0.15 M NaCl, oppure GBIS per test rapidi preliminari).

### Fase 2 — Benchmark Computazionale e Simulazione Pilot del Wild-Type (WT)
- [ ] **2.1** Eseguire un **benchmark di velocità** con **NAMD 3** (`namd3 +p<N_cores> ...`) per misurare la velocità in **`ns/day`** sul vostro hardware e annotare: numero di atomi del sistema, CPU/GPU usata, `ns/day`.
- [ ] **2.2** Definire e giustificare il **Computational Budget** complessivo ($\text{Sistemi} \times \text{Temperature} \times \text{Repliche} \times \text{Durata}$).
- [ ] **2.3** Eseguire la simulazione del **WT a 300 K** (controllo fisiologico) e ad **alta temperatura (es. 400 K)** (per indurre stress termico locale senza denaturare completamente il solenoide).
- [ ] **2.4** Analizzare la traiettoria WT con VMD:
  - Ricostruire le condizioni periodiche al contorno (`pbc wrap`) e allineare i frame sulla struttura iniziale (`frame 0`).
  - Calcolare **RMSD** del backbone e **RMSF** per residuo (`Cα`).
  - Calcolare il **$Q$-value globale** e il **$Q$-value per singolo residuo ($Q_i$)** ad alta temperatura per individuare i residui localmente instabili (es. $Q_i < 0.6$ nella finestra finale, come in Bekker et al.).

### Fase 3 — Ipotesi di Mutazione e Screening Energetico con FoldX
- [ ] **3.1** Ispezionare visivamente in VMD l'intorno 3D dei residui instabili/flessibili emersi dalla Fase 2 (cercando clash, cariche isolate, cavità non ottimali o H-bond/salt bridges mancanti tra repeat adiacenti).
- [ ] **3.2** Eseguire `RepairPDB` di **FoldX** su `2MYO.pdb` per ottenere la struttura ottimizzata `2MYO_Repair.pdb`.
- [ ] **3.3** Eseguire `PositionScan` di **FoldX** sulle posizioni candidate per esplorare le sostituzioni amminoacidiche e raccogliere i valori di $\Delta\Delta G$ ($\text{kcal/mol}$).
- [ ] **3.4** Selezionare le migliori mutazioni stabilizzanti ($\Delta\Delta G < 0$, idealmente $< -1.0\text{ kcal/mol}$) che siano strutturalmente sensate e lontane dai siti funzionali.
- [ ] **3.5** Costruire con `BuildModel` di **FoldX** i modelli 3D dei mutanti selezionati:
  - Almeno **2 mutanti singoli** (es. `Mut1`, `Mut2`) + **1 doppio mutante** (`Mut1 + Mut2`, seguendo l'esempio di Bekker et al. su `4idl`).

### Fase 4 — Simulazioni MD Comparative ("Matched Comparison")
- [ ] **4.1** Preparare tutti i sistemi (WT e varianti mutate) con **esattamente lo stesso protocollo** (stesso force field, stesso box/solvente, stessa concentrazione ionica).
- [ ] **4.2** Eseguire Minimizzazione + Equilibrazione + Produzione MD per **WT**, **Mutante 1**, **Mutante 2** e **Doppio Mutante** con **più repliche indipendenti** (cambiando il parametro `seed` delle velocità iniziali in NAMD, es. 3 repliche per sistema).

### Fase 5 — Analisi dei Risultati, Statistica e Presentazione
- [ ] **5.1** Calcolare il **$Q$-value** nel tempo e mediarlo sulla parte finale equilibrata della traiettoria (es. ultimi 30% della run) per tutte le repliche, riportando **Media $\pm$ Deviazione Standard** (come in Tabella II di Bekker et al.).
- [ ] **5.2** Opzionale ma consigliatissimo (punto bonus dal paper di Bekker et al.!): calcolare anche il $Q$-value diviso per gruppi chimico-fisici (**`all-all`**, **`hydrophilic-all`**, **`hydrophobic-small`**) e il numero di residui instabili residui in ciascun mutante.
- [ ] **5.3** Calcolare le osservabili strutturali di supporto: **RMSD**, **RMSF comparativo** (WT vs Mutanti sullo stesso grafico), **Hydrogen Bonds / Salt Bridges** formati dal residuo mutato, conservazione della struttura secondaria ($\alpha$-helicity).
- [ ] **5.4** Preparare le slide della presentazione (max 10 minuti) includendo il contributo specifico dei due membri del gruppo.

---

## 📂 4. Struttura della Repository `BML` e Organizzazione del Lavoro in Due

La repository è stata strutturata in cartelle numerate seguendo il flusso logico del progetto:

```text
BML/
├── README.md                  # Questo documento (Guida completa, Checklist e Roadmap)
├── .gitignore                 # Esclude file pesanti (.dcd, .coor, .vel, .log) da Git
├── 01_structures/             # Struttura PDB originale (2MYO), file puliti e PSF/PDB iniziali
├── 02_wt_baseline/            # Benchmark NAMD, simulazioni pilot WT (300K/400K) e individuazione siti instabili
├── 03_foldx_screening/        # RepairPDB, PositionScan, BuildModel e tabella dei ΔΔG predetti
├── 04_md_simulations/         # Configurazioni NAMD (.conf) per il confronto controllato WT vs Mutanti (multi-seed)
├── 05_analysis/               # Script TCL/Python, output Q-value, RMSD, RMSF, grafici e tabelle statistiche
├── 06_presentation/           # Slide, figure ad alta risoluzione da VMD e scaletta per il 16 Dicembre
└── scripts/                   # Script riutilizzabili (psfgen, calcolo Q-value, RMSF, plotting)
```

### 🤝 Come dividersi il lavoro in due (Oleg & Alessia) e usare Git senza conflitti
Alla fine della presentazione è richiesta una slide con la **"Brief description of each group member's contribution"**. Ecco una divisione del lavoro bilanciata che vi permette di lavorare in parallelo:

1. **Fase 1 & 3 (Insieme / Brainstorming condiviso):**
   - Analisi strutturale e funzionale di `2MYO`, scelta del protocollo MD dopo il benchmark, analisi dei residui instabili del WT e selezione razionale delle mutazioni da FoldX.
2. **Divisione Operativa durante Setup, Screening e MD (Fase 2, 3, 4, 5):**
   - **Membro 1 (es. Oleg):**
     - Preparazione della pipeline `psfgen` / solvatazione e benchmark computazionale iniziale (`ns/day`).
     - Esecuzione e analisi della simulazione WT di controllo (300 K) e del **Mutante Singolo 1** (tutte le repliche/seeds).
     - Sviluppo e automazione degli script TCL/Python per il calcolo del **$Q$-value** (globale, per gruppo e per residuo).
   - **Membro 2 (es. Alessia):**
     - Esecuzione dello screening sistematico con **FoldX** (`RepairPDB`, `PositionScan`, `BuildModel`) e verifica strutturale dei rotameri.
     - Esecuzione della simulazione WT ad alta temperatura (400 K), del **Mutante Singolo 2** e del **Doppio Mutante** (o divisione equa delle repliche sulla seconda macchina).
     - Analisi di **RMSD, RMSF, legami a idrogeno/ponti salini** e rendering grafico 3D in VMD dei contatti locali migliorati.

> [!IMPORTANT]
> **Regola Git per evitare problemi:** Non fate mai `git add` di file `.dcd`, `.coor` o `.vel` (sono già bloccati dal `.gitignore`). Condividete via Git gli script `.tcl`, i file `.conf` di NAMD, i `.pdb`/`.psf` di partenza e i file di testo dei risultati (`.dat`, `.txt`, `.csv`, `.fxout`).

---

## 🛠️ 5. Mappa Completa dei Tool e dei File del Corso

Tutti gli strumenti visti a lezione sono già installati sul vostro sistema o presenti nella cartella `/home/olly/UNI/BioMolEng2026`. Ecco la mappa completa con i riferimenti esatti:

| Tool / Risorsa | Percorso sul Sistema | A cosa serve nel Progetto |
| :--- | :--- | :--- |
| **NAMD 3.0.3** (Multicore) | `/usr/local/bin/namd3` | Motore di Dinamica Molecolare per minimizzazione, equilibrazione e produzione MD |
| **VMD 2.0.0** | `/usr/local/bin/vmd` | Visualizzazione 3D, preparazione topologia (`psfgen`, `solvate`, `autoionize`) e analisi TCL |
| **VMDextensions** (T. Giorgino) | `/usr/local/lib/vmd/plugins/VMDextensions.tcl` *(già caricato da `~/.vmdrc`)* | Libreria TCL che fornisce `prepareNativeContacts`, `measureNativeContacts`, `forFrames`, `qRMSD`, `helicity` |
| **Force Field CHARMM** | `/home/olly/UNI/BioMolEng2026/Tutorials/charmm36/` | File di topologia (`top_all36_prot.rtf`, `top_all27_prot_lipid.inp`) e parametri (`par_all36m_prot.prm`, `par_all27_prot_lipid.inp`) |
| **Tutorial Native Contacts ($Q$)** | `/home/olly/UNI/BioMolEng2026/Tutorials/nativecontacts/` | Contiene `native_contacts.tcl`, `psfgen.tcl` e `namd.conf` (esempio completo di calcolo di $Q$) |
| **Tutorial FoldX** | `/home/olly/UNI/BioMolEng2026/Tutorials/foldx_test/` | Esempi di `position_scan/`, `build_model/` e `stability/` con file `config.cfg` e `individual_list.txt` |
| **Binario FoldX** *(da recuperare da VM o sito CRG)* | *(Nella VM del corso era in `/home/biomolecular/Software/foldx`)* | Screening energetico $\Delta\Delta G$ (`RepairPDB`, `PositionScan`, `BuildModel`). *(Nota: per la sola costruzione 3D dei mutanti si può usare anche `mutate` in `psfgen` / VMD Mutator Plugin)* |
| **Script RMSF & RMSD per residuo** | `/home/olly/UNI/BioMolEng2026/Tutorials/rmsf.tcl` e `.../namd-tutorial-files/2-1-rmsd/residue_rmsd.tcl` | Script TCL per calcolare la fluttuazione quadratica media (RMSF) dei $C_\alpha$ e l'RMSD per residuo |
| **Configurazioni NAMD di riferimento** | `/home/olly/UNI/BioMolEng2026/Tutorials/namd-tutorial-files/1-3-box/ubq_wb_eq.conf` e `1-4-gbis/ubq_gbis_eq.conf` | Template NAMD per simulazioni in **Water Box esplicito (PME)** e in **Solvente Implicito (GBIS)** |
| **Benchmark NAMD** | `/home/olly/UNI/BioMolEng2026/Tutorials/dhfr/` | Esempio di benchmark per leggere i `ns/day` dai file di log di NAMD |

---

## 🗺️ 6. Roadmap Dettagliata Passo-Passo: Come Fare Ogni Fase

Di seguito trovate la spiegazione dettagliata di **cosa fare**, **perché farlo** e **quali comandi/tool usare** in ciascuna delle 5 fasi del workflow del progetto (Slide 15 di `Laboratory_Project_2026.pdf`).

---

### 📍 FASE 1: Analisi della Proteina Target e Preparazione della Struttura (`01_structures/`)

#### 1.1 Download e Ispezione Preliminare
Scaricate il file PDB di `2MYO` direttamente nella cartella `01_structures/`:
```bash
wget https://files.rcsb.org/download/2MYO.pdb -O 01_structures/2MYO.pdb
```
Apritelo in **VMD**:
```bash
vmd 01_structures/2MYO.pdb
```
**Cosa controllare in VMD (e documentare per le slide):**
1. **Rappresentazione grafica:** Impostate `Drawing Method: NewCartoon` e `Coloring Method: SecondaryStructure` (oppure `Index` per vedere la progressione N-term $\rightarrow$ C-term dei 4 Ankyrin Repeats).
2. **Completezza della sequenza:** Dal terminale TkConsole di VMD (`Extensions -> Tk Console`), grazie a `VMDextensions.tcl`, potete stampare la sequenza FASTA e controllare che non ci siano discontinuità (`/`):
   ```tcl
   set prot [atomselect top "protein"]
   getFasta $prot
   ```
   *(Verificherete che sono presenti tutti i 118 residui da `resid 1` a `resid 118` senza gap interni).*
3. **Stato delle Istidine (`HIS`):** In CHARMM le istidine neutre devono essere rinominate `HSD` (protonata su $N_\delta$), `HSE` (protonata su $N_\epsilon$) o `HSP` (carica $+1$, protonata su entrambi). Nei tutorial del corso (`Tutorials/charmm36/psfgen.tcl`) si usa l'alias standard `pdbalias residue HIS HSE`, ma controllate visivamente le istidine di `2MYO` in VMD per vedere se formano legami a idrogeno specifici.

#### 1.2 Generazione di `.psf` e `.pdb` con `psfgen`
Prendendo come riferimento `/home/olly/UNI/BioMolEng2026/Tutorials/nativecontacts/psfgen.tcl`, in VMD (o in modalità testuale `vmd -dispdev text -e script.tcl`) si estrae la sola proteina (senza gli idrogeni NMR vecchi per evitare conflitti di nomi atomici con CHARMM, usando `"chain A and protein and noh"`) e si genera la topologia:

```tcl
# Esempio concettuale basato su Tutorials/nativecontacts/psfgen.tcl
mol new 2MYO.pdb
set a [atomselect top "protein and noh"]
$a writepdb 2myo_clean.pdb

package require psfgen
topology /home/olly/UNI/BioMolEng2026/Tutorials/charmm36/top_all27_prot_lipid.inp
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
segment A {pdb 2myo_clean.pdb}
coordpdb 2myo_clean.pdb A
guesscoord
writepdb 2myo_wt.pdb
writepsf 2myo_wt.psf
```
*(Nota: potete usare `top_all27_prot_lipid.inp` + `par_all27_prot_lipid.inp` come negli script del corso, oppure il più recente `top_all36_prot.rtf` + `par_all36m_prot.prm` presente sempre in `Tutorials/charmm36/`).*

#### 1.3 Scelta del Modello di Solvente: Esplicito (Water Box + PME) vs Implicito (GBIS)
A lezione (`Lecture03`, Slide 35–41) e nei tutorial avete visto due approcci:
- **Opzione A — Solvente Esplicito TIP3P in Periodic Boundary Conditions (Raccomandata per fedeltà fisica, come in Bekker et al.):**
  - Si immerge `2myo_wt` in una scatola d'acqua TIP3P con un padding di almeno **10–12 Å** (Bekker et al. usano 15 Å, ma 10–12 Å riduce drasticamente il numero di atomi e accelera le run sul vostro PC!) e si aggiungono ioni $\text{Na}^+/\text{Cl}^-$ per neutralizzare il sistema e raggiungere una concentrazione fisiologica (**0.10 – 0.15 M NaCl**):
    ```tcl
    package require solvate
    solvate 2myo_wt.psf 2myo_wt.pdb -t 12 -o 2myo_wb

    package require autoionize
    autoionize -psf 2myo_wb.psf -pdb 2myo_wb.pdb -sc 0.15 -o 2myo_ionized
    ```
  - **Come calcolare le dimensioni del box per il file `.conf` di NAMD (`cellBasisVector1/2/3` e `cellOrigin`):**
    ```tcl
    set all [atomselect top all]
    measure minmax $all
    measure center $all
    ```
    La differenza `max - min` lungo X, Y, Z vi dà i tre `cellBasisVector`, mentre `measure center $all` vi dà `cellOrigin`!
- **Opzione B — Solvente Implicito GBIS (`gbis on`, come in `Tutorials/nativecontacts/namd.conf` e `1-4-gbis/ubq_gbis_eq.conf`):**
  - Molto più veloce perché simula solo i ~1850 atomi della proteina senza molecole d'acqua esplicite, ma (come spiegato in `Lecture03`, Slide 40) approssima l'effetto idrofobico, non modella i ponti idrogeno specifici col solvente e manca della viscosità dell'acqua (quindi ad alte temperature la proteina si deforma molto più velocemente!).
  - **Consiglio strategico:** Fate subito il benchmark (Fase 2.1) con il box d'acqua esplicito. Con una proteina piccola di 118 aa e NAMD 3 multicore, il solvente esplicito (o un padding compatto di 10–12 Å) è solitamente gestibilissimo e scientificamente molto più solido per valutare mutazioni di residui idrofilici/idrofobici!

---

### 📍 FASE 2: Benchmark, Computational Budget e Caratterizzazione del WT (`02_wt_baseline/`)

#### 2.1 Benchmark di Velocità (`ns/day`) e Pianificazione del Budget
Prima di lanciare simulazioni lunghe, fate una breve run di test (es. 5.000 – 10.000 step = 10–20 ps) con **NAMD 3** usando tutti i core fisici/logici della vostra CPU:
```bash
namd3 +p8 test_benchmark.conf > test_benchmark.log
```
Aprite `test_benchmark.log` e cercate le righe `Benchmark time:` (oppure `TIMING:`): NAMD stampa direttamente i **`days/ns`** (da cui $\text{ns/day} = \frac{1}{\text{days/ns}}$).
- **Calcolate e documentate il vostro Computational Budget (Slide 21 di `Laboratory_Project_2026.pdf`):**
  $$\text{Total MD Time} = N_{\text{systems}} \times N_{\text{temperatures}} \times N_{\text{replicas (seeds)}} \times t_{\text{trajectory}}$$
  *Esempio pratico:* Se misurate ad esempio $20\text{ ns/day}$ su un PC (e lavorate su 2 PC in parallelo = $40\text{ ns/day}$ complessivi), un piano come **4 sistemi** (WT, Mut1, Mut2, Doppio Mutante) $\times$ **1 temperatura produttiva di confronto** (es. 400 K, più il controllo WT a 300 K) $\times$ **3 repliche (seeds)** $\times$ **15–20 ns** è perfettamente realizzabile in pochi giorni di calcolo ed è statisticamente rigoroso!

#### 2.2 Simulazione Pilot del Wild-Type (300 K vs 400 K)
Nel paper di Bekker et al.:
- A **300 K** la struttura resta molto vicina a quella nativa ($Q \approx 0.95 - 0.97$) e discrimina poco la stabilità ($r = 0.51$).
- A **400 K** le fluttuazioni termiche sono accelerate ma il fold globale non crolla: i residui più deboli perdono i contatti nativi e il $Q$-value mostra una forte correlazione con la $T_m$ sperimentale ($r = 0.79$, che sale a $r = 0.84$ guardando i contatti `hydrophilic-all`).
- A **500 K** la perturbazione è troppo violenta e tutte le proteine perdono rapidamente il fold ($r = 0.08$).

> [!TIP]
> **Verifica sulla Myotrophin (`2MYO`):** Essendo `2MYO` una proteina a 4 Ankyrin Repeats (puramente $\alpha$-elica, senza ponti disolfuro, a differenza del $\beta$-sandwich degli anticorpi di Bekker et al. che hanno un ponte disolfuro interno), verificate nella pilot run del WT come si comporta a **400 K** (ed eventualmente confrontate con **350 K – 370 K** se a 400 K dovesse srotolarsi troppo in fretta, oppure confermate **400 K** se mantiene il core $\alpha$-elicoidale perdendo solo contatti locali!). Come sottolinea la Slide 12 del progetto: *"A useful temperature provides discrimination without making every system rapidly lose its fold."*

#### 2.3 Come individuare le regioni e i residui instabili del WT
Dalla traiettoria pilot del WT ad alta temperatura, usate due metriche complementari:

1. **RMSF (Root Mean Square Fluctuation) per residuo (Slide 22):**
   - **Attenzione obbligatoria prima del calcolo:** In VMD dovete prima ricomporre le molecole attraverso i bordi periodici e **allineare tutti i frame sulla struttura di riferimento** (altrimenti calcolate la traslazione/rotazione della proteina nella scatola invece della flessibilità interna!):
     ```tcl
     # 1. Ricostruzione PBC (se in water box)
     package require pbctools
     pbc wrap -compound residue -center com -centersel "protein" -all

     # 2. Allineamento di tutti i frame sul frame 0 (backbone o CA)
     set ref [atomselect top "protein and name CA" frame 0]
     set sel [atomselect top "protein and name CA"]
     set all [atomselect top all]
     for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
         $sel frame $i
         $all frame $i
         $all move [measure fit $sel $ref]
     }

     # 3. Calcolo RMSF (vedi /home/olly/UNI/BioMolEng2026/Tutorials/rmsf.tcl)
     set outfile [open rmsf_wt.dat w]
     set rmsf [measure rmsf $sel]
     for {set i 0} {$i < [$sel num]} {incr i} {
         puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
     }
     close $outfile
     ```
   - I **picchi di RMSF** indicano le regioni con maggiore fluttuazione posizionale. Ricordate però la nota della Slide 22: *"flexibility alone does not establish local instability or identify a beneficial mutation"* (es. N-term, C-term e loop funzionali sono naturalmente flessibili!).

2. **$Q$-value per singolo residuo ($Q_i$) secondo Bekker et al.:**
   - Per ogni residuo $i$ (`resid 1` ... `118`), possiamo usare `prepareNativeContacts` e `measureNativeContacts` di `VMDextensions.tcl` passando come `sel1` il residuo `resid $i and noh` e come `sel2` il resto della proteina (`protein and noh and not resid $i`), mediando sugli ultimi ns della simulazione:
     ```tcl
     # Calcolo del Q-value per singolo residuo i (approccio Bekker et al.)
     set prot [atomselect top "protein and noh"]
     for {set r 1} {$r <= 118} {incr r} {
         set res_sel [atomselect top "protein and noh and resid $r" frame 0]
         set env_sel [atomselect top "protein and noh and not resid $r" frame 0]
         set ref_r [prepareNativeContacts 4.5 $res_sel $env_sel]
         set nnc_r [llength $ref_r]
         # Se nnc_r >= 30 (soglia di Bekker et al.), misurare measureNativeContacts sugli ultimi frame
         # I residui con <Q_r> < 0.6 sono classificati come "unstable residues"!
     }
     ```
   - Questo vi darà una **tabella identica alla Table S1 di Bekker et al.**, identificando con precisione quali residui (al di fuori di N/C-term e dei loop di legame con CP) perdono i propri contatti nativi ad alta temperatura!

---

### 📍 FASE 3: Progettazione delle Mutazioni e Screening con FoldX (`03_foldx_screening/`)

Una volta individuati i siti candidati (residui con basso $Q_i$ o alto RMSF non funzionale, o contatti sub-ottimali tra gli Ankyrin Repeats), usate **FoldX** per progettare e valutare energeticamente le mutazioni (Slide 16–19 e `/home/olly/UNI/BioMolEng2026/Tutorials/foldx_test/`).

> [!NOTE]
> **Criteri biofisici per proporre mutazioni negli Ankyrin Repeats di `2MYO`:**
> - **Consenso degli Ankyrin Repeats:** Le proteine ad Ankyrin Repeat hanno una forte regolarità di sequenza tra i vari repeat (motivi conservati come il tetrapeptide `TPLH`, residui idrofobici all'interfaccia tra le eliche antiparallele, e interazioni elettrostatiche/H-bond tra repeat adiacenti). Confrontare i 4 repeat di `2MYO` con il consenso degli ankyrin repeat o analizzare dove un repeat devia dal consenso è un'ottima giustificazione scientifica per proporre una mutazione!
> - **Interazioni di residui idrofilici (Bekker et al.):** Bekker et al. dimostrano che la stabilità termica correla fortemente con il mantenimento dei contatti dei **residui idrofilici** (`Asp, Glu, Gln, Asn, Arg, Lys, His`) e il corretto posizionamento di residui idrofobici/idrofilici (es. nel loro caso studio su `4idl`, sostituiscono `R71I` in un loop per migliorare il packing idrofobico e `N27D` per formare interazioni elettrostatiche più forti).

####Il Workflow FoldX in 3 Passi:

1. **Passo 1 — `RepairPDB` (Obbligatorio prima di ogni analisi FoldX):**
   Rimuove i clash sterici locali e ottimizza i rotameri delle catene laterali del PDB di partenza:
   ```bash
   foldx --command=RepairPDB --pdb=2MYO.pdb
   ```
   *(Genera `2MYO_Repair.pdb`, che sarà l'input per tutti i comandi successivi).*

2. **Passo 2 — `PositionScan` (Esplorazione delle sostituzioni nei siti scelti):**
   Crea un file di configurazione `config_ps.cfg` (vedi `Tutorials/foldx_test/position_scan/config.cfg`) oppure lancialo da riga di comando.
   Il formato della posizione in FoldX è: **`<Residuo 1-lettera><Catena><NumeroResiduo><CodiceMutazione>`**, dove il codice mutazione può essere:
   - `a` = tutti i 20 amminoacidi naturali
   - `h` = residui idrofobici (`{GLY, ALA, LEU, VAL, ILE, THR, SER, CYS, MET, LYS, TRP, TYR, PHE, HIS}`)
   - `c` = residui carichi (`{ARG, LYS, GLU, ASP, HIS}`)
   - `p` = residui polari (`{ARG, LYS, GLU, ASP, HIS, TRP, TYR, THR, SER, GLN, ASN}`)
   - oppure la singola lettera dell'amminoacido target.
   ```ini
   command=PositionScan
   pdb=2MYO_Repair.pdb
   positions=RA40a,KA65a
   ```
   Guardate il file di output `PS_2MYO_Repair_scanning_output.txt` e cercate le mutazioni con **$\Delta\Delta G < 0\text{ kcal/mol}$** (più è negativo, maggiore è la stabilizzazione predetta da FoldX).

3. **Passo 3 — `BuildModel` (Costruzione dei mutanti singoli e combinati/doppi):**
   Per i candidati scelti (es. 2 mutazioni singole e la loro combinazione doppia), create il file obbligatorio **`individual_list.txt`** (attenzione: ogni riga rappresenta un mutante e **deve terminare con `;`**!):
   ```text
   RA40E;
   KA65I;
   RA40E,KA65I;
   ```
   *(Dove `R` = residuo WT, `A` = catena A, `40` = numero residuo, `E` = nuovo residuo).*
   Lanciate `BuildModel` (vedi `Tutorials/foldx_test/build_model/config.cfg`):
   ```ini
   command=BuildModel
   pdb=2MYO_Repair.pdb
   mutant-file=individual_list.txt
   numberOfRuns=3
   ```
   - **Output di `BuildModel`:**
     - `Dif_2MYO_Repair.fxout` e `Average_2MYO_Repair.fxout`: contengono il **$\Delta\Delta G$ totale** (`total energy`) e la scomposizione in tutti i termini energetici (*Backbone Hbond, Sidechain Hbond, Van der Waals, Electrostatics, Solvation Polar, Solvation Hydrophobic, Van der Waals clashes, entropy...*) $\rightarrow$ **Ottimo da inserire in una tabella nelle slide per spiegare *perché* la mutazione stabilizza la proteina!**
     - `2MYO_Repair_1.pdb`, `2MYO_Repair_2.pdb`, `2MYO_Repair_3.pdb`: i PDB dei vostri mutanti pronti per essere passati a `psfgen` e simulati in NAMD!
     - `WT_2MYO_Repair_1.pdb`, ecc.: la struttura di riferimento WT ottimizzata sullo stesso intorno locale.

---

### 📍 FASE 4: Simulazioni MD Comparative ("Matched Comparison", `04_md_simulations/`)

Questa è la fase centrale di validazione dinamica (Slide 15 & 21):
1. **Protocollo "Matched" (Controllo rigoroso delle variabili):**
   - Costruite i sistemi per **WT**, **Mutante 1**, **Mutante 2** e **Doppio Mutante** usando **esattamente lo stesso script `psfgen` + `solvate` + `autoionize`** e **lo stesso file di configurazione `.conf` di NAMD** (stessa temperatura, stesso cutoff, stesso timestep di `2.0 fs` con `rigidBonds all`, stesso termostato `langevin on`).
2. **Protocollo di Equilibrazione e Produzione (vedi `Lecture04`, Slide 26–29 e `Tutorials/namd-tutorial-files/1-3-box/ubq_wb_eq.conf`):**
   - **Minimizzazione energetica:** `minimize 1000` (o `100`–`500` step) per rimuovere eventuali clash di idrogeni/acqua.
   - **Equilibrazione:** Riscaldamento/rilassamento del solvente in NPT a 1 atm (`langevinPiston on`, `langevinPistonTarget 1.01325`) per stabilizzare la densità del box d'acqua (~1.0 g/cm³).
   - **Produzione ad alta temperatura (es. 400 K):** Run non vincolata per la durata stabilita nel vostro computational budget, salvando le coordinate nel file `.dcd` ad intervalli regolari (es. ogni 5–10 ps, ossia `dcdfreq 2500` o `5000` con timestep di 2 fs).
3. **Repliche indipendenti tramite `seed` delle velocità iniziali:**
   - Per avere barre d'errore (media $\pm$ deviazione standard) come nel paper di Bekker et al., lanciate **almeno 3 repliche indipendenti** per ciascun sistema cambiando in NAMD il parametro:
     ```tcl
     seed 12345   ;# Cambiare numero intero positivo per Replica 1, Replica 2, Replica 3
     ```
     Questo inizializza le velocità atomiche dalla distribuzione di Maxwell-Boltzmann con numeri casuali diversi, generando traiettorie indipendenti.

---

### 📍 FASE 5: Analisi del $Q$-value e delle Osservabili di Supporto (`05_analysis/`)

#### 5.1 Calcolo del $Q$-value (Fraction of Native Contacts) con `VMDextensions.tcl`
Il $Q$-value misura la frazione di contatti atomici presenti nello stato nativo (frame 0 / struttura iniziale) che vengono mantenuti nel corso della simulazione a temperatura elevata:
$$Q(t) = \frac{N_{\text{native contacts}}(t)}{N_{\text{native contacts}}(0)} \in [0, 1]$$
Un valore di **$Q$ più alto** negli ultimi nanosecondi di simulazione indica una **maggiore ritenzione della struttura nativa** e quindi una **maggiore stabilità termica**.

Prendendo come base `/home/olly/UNI/BioMolEng2026/Tutorials/nativecontacts/native_contacts.tcl`, potete calcolare:
1. **$Q$-value globale ($C_\alpha$ con cutoff 7.0 Å come nel tutorial del corso, oppure atomi pesanti `noh` con cutoff 4.5 Å):**
   ```tcl
   # Esempio basato su Tutorials/nativecontacts/native_contacts.tcl
   set protein [atomselect top "protein and name CA"]
   $protein frame 0
   set ref [prepareNativeContacts 7.0 $protein]
   set nnc [llength $ref]
   set out [open q_global.dat w]

   forFrames fno $protein {
       set nc [measureNativeContacts $ref 7.0 $protein]
       set qnc [expr {1.0 * $nc / $nnc}]
       puts $out "$fno $qnc"
   }
   close $out
   ```
2. **$Q$-value per classi chimico-fisiche di residui (come in Figura 2 e Tabella II di Bekker et al.!):**
   Nel paper di riferimento, gli autori classificano gli amminoacidi in 3 gruppi:
   - **Hydrophilic:** `resname ASP GLU GLN ASN ARG LYS HIS HSE HSD HSP`
   - **Hydrophobic:** `resname PHE TYR TRP LEU VAL ILE MET CYS PRO`
   - **Small:** `resname GLY ALA SER THR`
   Passando due selezioni a `prepareNativeContacts` e `measureNativeContacts`, potete calcolare non solo `Q: all-all`, ma anche **`Q: hydrophilic-all`** (che nel paper ha la correlazione più alta in assoluto con la $T_m$, $r = 0.84$!) e **`Q: hydrophobic-small`**:
   ```tcl
   set hydro [atomselect top "protein and noh and resname ASP GLU GLN ASN ARG LYS HIS HSE HSD HSP"]
   set allp  [atomselect top "protein and noh"]
   $hydro frame 0; $allp frame 0
   set ref_ha [prepareNativeContacts 4.5 $hydro $allp]
   # Poi nel ciclo forFrames aggiornate sia $hydro che $allp ($allp frame $fno) e misurate:
   # measureNativeContacts $ref_ha 4.5 $hydro $allp
   ```

#### 5.2 Media Statistica e Confronto tra WT e Varianti
- Per ogni replica (seed), calcolate la media temporale di $Q(t)$ e dell'RMSD sulla **finestra finale della traiettoria** (es. gli ultimi 30 ns su 100 ns in Bekker et al., oppure l'ultimo 30% della vostra durata di produzione, quando il decadimento iniziale si è assestato).
- Calcolate poi **Media $\pm$ Deviazione Standard** tra le repliche e costruite una tabella riassuntiva identica alla **Table II di Bekker et al.**:

| Sistema Computazionale | $\Delta\Delta G_{\text{FoldX}}$ (kcal/mol) | RMSD finale (Å) | $Q$: `all-all` | $Q$: `hydrophilic-all` | $Q$: `hydrophobic-small` |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **2MYO Wild-Type (WT)** | — | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ |
| **Mutante 1** | ... | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ |
| **Mutante 2** | ... | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ |
| **Doppio Mutante (1+2)** | ... | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ | $\mu \pm \sigma$ |

#### 5.3 Osservabili di Supporto per spiegare il meccanismo fisico
Per ottenere il massimo punteggio nel **Learning Objective 4 (Scientific Interpretation)**, affiancate al $Q$-value:
1. **Confronto RMSF per residuo (WT vs Mutanti):** Sovrapponete la curva RMSF del WT e del Mutante (proprio come nella Slide 34 di `Lecture04`!) per mostrare visivamente l'abbassamento del picco di flessibilità nella regione mutata.
2. **Occupazione di Hydrogen Bonds e Salt Bridges (`VMD -> Extensions -> Analysis -> Hydrogen Bonds / Salt Bridges`):** Calcolate la `% di occupancy` del nuovo legame a idrogeno o ponte salino introdotto dalla vostra mutazione durante la traiettoria MD.
3. **Frazione di $\alpha$-elica (`helicity` da `VMDextensions.tcl` o VMD Timeline):** Essendo `2MYO` composta da 4 Ankyrin Repeats ad $\alpha$-elica, potete monitorare in una riga con `helicity [atomselect top "name CA"]` quanto bene viene preservato il contenuto elicoidale ad alta temperatura nel mutante rispetto al WT!

---

## 🎤 7. Struttura Consigliata per la Presentazione Finale (16 Dicembre — Max 10 Minuti)

Con 10 minuti a disposizione in due (~5 minuti a testa), vi conviene preparare circa **8–10 slide** molto visive e dritte al punto:

1. **Slide 1 — Titolo & Team #10:** Titolo del progetto, nomi (Oleg Nedina, Alessia Rigoni), proteina assegnata (`2MYO` - Myotrophin).
2. **Slide 2 — Introduzione alla Myotrophin (`2MYO`):** Architettura a 4 Ankyrin Repeats (struttura NMR, 118 aa), folding cooperativo a 2 stati, funzione biologica (inibizione dell'Actin-Capping Protein e ruolo nell'ipertrofia cardiaca) e **vincoli funzionali** (quali loop/residui non devono essere toccati).
3. **Slide 3 — Metodologia & Computational Budget:** Schema del workflow (Baseline WT $\rightarrow$ FoldX Screening $\rightarrow$ Matched MD $\rightarrow$ Analisi $Q$-value), benchmark `ns/day` misurato, protocollo MD scelto (temperatura, solvente, numero di repliche/seeds e durata giustificati dal budget).
4. **Slide 4 — Caratterizzazione Baseline del WT (300 K vs Alta Temperatura):** Grafici di RMSD, RMSF e $Q_i$ per residuo del WT ad alta temperatura; identificazione strutturale 3D dei residui localmente instabili nei 4 Ankyrin Repeats.
5. **Slide 5 — Razionale delle Mutazioni e Screening FoldX:** Ipotesi chimico-fisiche sulle posizioni scelte, risultati di `PositionScan` e `BuildModel` ($\Delta\Delta G$ e contributi energetici di H-bond/VdW/solvatazione), immagini zoomate in VMD del sito prima (WT) e dopo la mutazione.
6. **Slide 6 — Risultati MD Comparativi ($Q$-value nel tempo e Tabella Statistica):** Curve di $Q(t)$ e Tabella riassuntiva (Media $\pm$ SD su tutte le repliche per WT, mutanti singoli e doppio mutante, inclusi `all-all` e `hydrophilic-all`).
7. **Slide 7 — Evidenze Strutturali di Supporto (Meccanismo di Stabilizzazione):** Confronto RMSF (riduzione del picco locale nel mutante), occupancy di H-bond/Salt Bridges e conservazione delle $\alpha$-eliche.
8. **Slide 8 — Conclusioni, Limiti del Modello e Contributi del Team:** Discussione critica dell'incertezza statistica tra le repliche, come questo screening computazionale guida la validazione sperimentale *in vitro* (CD spectroscopy / misura di $T_m$), e **riepilogo del contributo di Oleg e Alessia**.

---

## 📚 8. Riferimenti Bibliografici Utili per il Progetto

1. **Studio di Riferimento del Corso:**  
   Bekker GJ, Ma B, Kamiya N (2019). *Thermal stability of single-domain antibodies estimated by molecular dynamics simulations.* **Protein Science**, 28(2):429–438. [doi:10.1002/pro.3546](https://doi.org/10.1002/pro.3546)
2. **Struttura NMR della Myotrophin (`2MYO` / `1MYO`):**  
   Yang Y, Nanduri S, Sen S, Qin J (1998). *The structural basis of ankyrin-like repeat function as revealed by the solution structure of myotrophin.* **Structure**, 6(5):619–626. [RCSB: 2MYO](https://www.rcsb.org/structure/2MYO)
3. **Complesso Funzionale Myotrophin (V-1) / Actin-Capping Protein (`3AAA`):**  
   Takeda S, Minakata S, Koike R, Kawahata I, Narita A, Kitazawa M, Ota M, Yamakuni T, Maéda Y, Nitanai Y (2010). *Two distinct mechanisms for actin capping protein regulation—steric and allosteric inhibition.* **PLoS Biology**, 8(7):e1000416.
4. **Libreria VMDextensions ($Q$-value / Native Contacts):**  
   Giorgino T (2019). *Analysis Libraries for Molecular Trajectories: A Cross-Language Synopsis.* In: **Biomolecular Simulations: Methods and Protocols** (Methods in Molecular Biology), Springer, pp. 503–527. [Documentazione Native Contacts](http://tonigi.github.io/vmd_extensions/group__native.html)
5. **FoldX Suite:**  
   Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L (2005). *The FoldX web server: an online force field.* **Nucleic Acids Research**, 33:W382–W388. [FoldX Documentation](http://foldxsuite.crg.eu/documentation)
