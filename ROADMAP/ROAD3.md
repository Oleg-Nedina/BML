# 🗺️ ROAD3.md — Fase 3: Multi-Tool Consensus Screening & Modellazione dei Mutanti (`03_mutations/`)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Target:** Myotrophin Wild-Type (`2MYO.pdb`, 118 aa, `Chain A`)  
**Cartella di lavoro di questa fase:** `03_mutations/`  
**Regole di Riferimento:** Tutte le operazioni seguono tassativamente quanto stabilito in [`AGENTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/AGENTS.md) (*Zero Ipotesi, Citazione Obbligatoria della Documentazione, Esecuzione da parte dell'Utente, Grafici in formato vettoriale `.pdf`, Note Foto per le Slide*).

---

## 📚 1. Fonti e Documentazione di Riferimento per la Fase 3

Ogni comando, parametro e criterio decisionale di questa fase deriva da:

1. **Linee Guida Ufficiali del Progetto:**  
   [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 14–16** (*Workflow Step 3: FoldX Mutation Screening — RepairPDB, PositionScan, BuildModel*), **Slide 17–20** (*Selection of 2 single mutations + 1 double mutant, avoiding active/binding sites*).
2. **Slide delle Lezioni del Corso (Prof. A. Gautieri):**
   - [`Lecture05 - Protein Engineering and Design.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture05%20-%20Protein%20Engineering%20and%20Design.pdf) — **Slide 20–32** (Approcci evolutivi *Back-to-Consensus* e allineamento dei moduli), **Slide 33–47** (Funzione di energia empirica di **FoldX**: $\Delta G = W_{\text{vdw}}\Delta G_{\text{vdw}} + W_{\text{solvH}}\Delta G_{\text{solvH}} + \dots$, `RepairPDB`, `PositionScan`, `BuildModel` e interpretazione di $\Delta\Delta G < -0.5 / -1.0\text{ kcal/mol}$), **Slide 48–59** (**Rosetta** e *Rotamer Packing*), **Slide 60–70** (Metodi AI Graph-based: **ProteinMPNN / ThermoMPNN** e approcci ibridi tipo **FireProt**).
3. **File dei Tutorial Ufficiali del Corso (`Tutorials/foldx_test/`):**
   - [`Tutorials/foldx_test/position_scan/config.cfg`](file:///home/olly/UNI/BioMolEng2026/Tutorials/foldx_test/position_scan/config.cfg) (`command=PositionScan`, sintassi `positions=<WT><Chain><ResNum><Target>`, dove `d` scansiona tutti i 20 amminoacidi).
   - [`Tutorials/foldx_test/build_model/config.cfg`](file:///home/olly/UNI/BioMolEng2026/Tutorials/foldx_test/build_model/config.cfg) (`command=BuildModel`, `mutant-file=individual_list.txt` con punto e virgola `;` a fine riga).
4. **Letteratura Scientifica sulle Ankyrin Repeat Proteins & Myotrophin:**
   - **Mosavi LK, Minor DL, Peng ZY (2002).** *Consensus-derived structural determinants of the ankyrin repeat motif.* **PNAS**, 99(25):16029–16034 (Sequenza consenso di 33 residui degli Ankyrin Repeats termostabili).
   - **Lowe AR, Itzhaki LS (2007).** *Rational redesign of the folding pathway of a modular protein (Myotrophin `2MYO`).* **PNAS**, 104(8):2679–2684.
   - **Bekker GJ, Ma B, Kamiya N (2019).** *Thermal stability of single-domain antibodies estimated by molecular dynamics simulations.* **Protein Science**, 28(2):429–438 (Incrocio tra residui dinamicamente instabili $\langle Q_i \rangle < 0.6$ a 400 K e design razionale).

---

## 🎯 2. L'Architettura della Fase 3 ("Orthogonal Consensus Funnel" spiegato a soldoni)

Perché non ci limitiamo a scegliere 2 mutazioni a caso con FoldX, ma usiamo un **"Imbuto a 4 Filtri Ortogonali"**?
- **Il limite di un singolo tool statico (es. solo FoldX):** FoldX calcola la differenza di energia libera ($\Delta\Delta G = \Delta G_{\text{Mut}} - \Delta G_{\text{WT}}$) tenendo lo scheletro della proteina (*backbone*) quasi fermo. A volte FoldX dà un voto ottimo ($\Delta\Delta G < -1.5\text{ kcal/mol}$) a una mutazione che infila un amminoacido enorme dentro un loop, ma appena la metti in acqua calda (Dinamica Molecolare a 400 K) quella mutazione fa irrigidire o spaccare il loop!
- **La nostra soluzione da ricerca (i 4 Pilastri che possiamo già calcolare OGGI mentre aspettiamo la MD notturna!):**

| Pilastro | Cosa fa (in termini di Ingegneria / HPC) | Quando lo eseguiamo? | Output generato |
| :--- | :--- | :--- | :--- |
| **Filtro 0: Functional Masking** | **Protezione dei pin di I/O (`3AAA`):** vieta qualsiasi modifica sui **24 residui** che toccano l'Actin-Capping Protein ($< 4.5\text{ Å}$). Restano **94 residui** modificabili. | ✅ **Già completato in Fase 1** | Lista dei 24 residui mascherati |
| **Pilastro 1: Riparazione & Scansione Sistematica FoldX 5** | 1) `RepairPDB` rilassa le catene laterali del PDB NMR `2MYO.pdb`.<br>2) `PositionScan` prova **tutti i 20 amminoacidi** sulle posizioni chiave (o su tutti i 94 residui dello scaffold!), calcolando la matrice completa di $\Delta\Delta G_{\text{FoldX}}$. | ⚡ **SUBITO (prima della MD notturna!)** | `03_mutations/2MYO_Repair.pdb` + classifica completa $\Delta\Delta G_{\text{FoldX}}$ |
| **Pilastro 2: Consenso Evolutivo degli Ankyrin Repeats (`ANK1–4`)** | Confronta i 4 mattoncini di `2MYO` contro il "mattoncino Ankyrin perfetto" di **Mosavi et al. (PNAS 2002)** + analizza gli angoli torsionali anomali del `REMARK 500` di `2MYO.pdb`. | ⚡ **SUBITO (0.1 secondi in Python!)** | `03_mutations/ankyrin_consensus_report.txt` + grafico `.pdf` |
| **Pilastro 3: Consenso Cloud Multi-Tool (FireProt / ThermoMPNN / DynaMut2)** | Invia `2MYO_Repair.pdb` ai server gratuiti **ThermoMPNN** (Graph Neural Network della `Lecture05`), **FireProt 2.0** (FoldX + Rosetta + Evoluzione) e **DynaMut2** (Normal Mode Analysis). | ⚡ **SUBITO (lavorano in Cloud mentre dormi!)** | Tabella $\Delta\Delta G_{\text{ThermoMPNN}}$, $\Delta\Delta G_{\text{Rosetta/FireProt}}$, $\Delta\Delta G_{\text{DynaMut2}}$ |
| **Pilastro 4: Instabilità Dinamica MD ($\langle Q_i \rangle_{400\text{K}} < 0.6$)** | La simulazione MD ad alta temperatura ($400\text{ K}$) ci dice quali residui **perdono davvero i contatti nativi** sotto stress termico. | 🌙 **Stanotte (Run MD) $\rightarrow$ Domani mattina incrociamo!** | `02_wt_baseline/unstable_residues_report.txt` |

> **💡 Il grande vantaggio di fare i Pilastri 1, 2 e 3 ADESSO:**  
> Domani mattina, appena finisce la simulazione MD del Wild-Type, avremo **già pronte tutte le tabelle energetiche di FoldX, del Consenso Ankyrin e dei server Cloud**. Ci basterà fare l'intersezione (in 2 minuti!) con i residui instabili emersi dalla MD ($\langle Q_i \rangle < 0.6$) e lanciare subito `BuildModel` per generare i 3 mutanti (`Mut1`, `Mut2`, `Double`)!

---

## 🛠️ PASSO 3.1 — Analisi Automatica del Consenso Ankyrin Repeat & Struttura NMR (`scripts/04_ankyrin_consensus_and_targets.py`)

Prima ancora di lanciare FoldX, facciamo girare uno script Python velocissimo (`scripts/04_ankyrin_consensus_and_targets.py`) che incrocia **3 informazioni strutturali già presenti nei nostri file**:
1. **Allineamento posizione per posizione (1..33)** dei 4 Ankyrin Repeats di `2MYO` (`ANK1`, `ANK2`, `ANK3`, `ANK4`) contro la sequenza consenso universale degli Ankyrin Repeats termostabili (**Mosavi et al., *PNAS* 2002**):
   $$\texttt{1:G 2:x 3:T 4:P 5:L 6:H 7:L 8:A 9:A 10:x 11:x 12:G 13:H 14:L 15:E 16:I 17:V 18:E 19:V 20:L 21:L 22:K 23:x 24:G 25:A 26:D 27:V 28:N 29:A 30:x 31:D 32:x 33:x}$$
2. **Estrazione delle anomalie geometriche dal `REMARK 500` di `01_structures/raw/2MYO.pdb`** (righe 104–160: i 4 contatti sterici ravvicinati e i 35 residui con angoli di Ramachandran $(\phi, \psi)$ fuori dalle regioni canoniche nell'NMR).
3. **Filtro automatico della Functional Mask (24 residui):** scarta automaticamente i 24 residui che legano l'Actin-Capping Protein in `3AAA` e genera la lista esatta delle posizioni non-mascherate per FoldX!

### 💻 Comando 3.1 — Eseguire l'analisi del Consenso Ankyrin e generare il grafico `.pdf`
Dalla cartella `BML/`:
```bash
python3 scripts/04_ankyrin_consensus_and_targets.py
```
- **Output numerico:** `03_mutations/ankyrin_consensus_report.txt` (tabella completa di tutti i 4 repeat, le deviazioni dal consenso come `L73 -> H`, `K94 -> V/I`, `P64 -> K/N`, `L105 -> F` Ratto/Uomo, ecc.).
- **Output grafico vettoriale:** `06_presentation/fig6_ankyrin_consensus_alignment.pdf` (mappa visuale a colori dei 4 repeat allineati con il consenso e con la Functional Mask, pronta per le Slide!).

---

## 🔧 PASSO 3.2 — Riparazione della Struttura WT con FoldX 5 (`RepairPDB`)

Come spiegato nella `Lecture05` (Slide 37) e nel tutorial `Tutorials/foldx_test/`, **non si lancia mai uno studio di mutazione direttamente su un PDB grezzo del Protein Data Bank**, specialmente se è una struttura NMR come `2MYO.pdb`.
Il comando **`RepairPDB`** di FoldX:
- Tiene fisso il backbone e **ottimizza i rotameri delle catene laterali** che hanno piccoli scontri di Van der Waals o legami a idrogeno non ottimali;
- Esporta la struttura rilassata **`03_mutations/2MYO_Repair.pdb`**, che sarà l'**unico input valido** per `PositionScan` e `BuildModel`.

### 💻 Comando 3.2 — Eseguire `RepairPDB` con FoldX 5 (~60 secondi)
Abbiamo preparato il file di configurazione `03_mutations/repair.cfg`. Lancialo dalla cartella `BML/`:
```bash
# Copia 2MYO.pdb in 03_mutations/ ed esegue RepairPDB
cp 01_structures/raw/2MYO.pdb 03_mutations/2MYO.pdb
foldx -f 03_mutations/repair.cfg > 03_mutations/repair.log
```
Quando termina (circa 1 minuto), verifica l'energia totale prima e dopo la riparazione:
```bash
grep -E "Total|Backbone|Sidechain|Van der Waals" 03_mutations/2MYO_Repair.fxout 03_mutations/repair.log | tail -n 20
```
- **Cosa aspettarsi:** Vedrai che `03_mutations/2MYO_Repair.pdb` è stato creato e che l'energia libera totale $\Delta G$ di FoldX è scesa di diversi $\text{kcal/mol}$ (eliminando gli scontri sterici dell'NMR!).

---

## 🔍 PASSO 3.3 — Scansione Sistematica con FoldX `PositionScan` sulle Posizioni Chiave

Come funziona **`PositionScan`** (da `Tutorials/foldx_test/position_scan/config.cfg`)?
- La sintassi per ogni posizione è: **`<AA_WT><Catena><NumeroResiduo><Target>`**
  - Esempio: **`LA73d`** significa: *"Prendi la **L**eucina (**L**) della catena **A** al residuo **73** e mútala in tutti i **20 amminoacidi standard** (**d**)"*.
- FoldX costruisce tutte le 20 varianti per ciascuna posizione, ottimizza i residui vicini e scrive nel file **`PS_2MYO_Repair_scanning_output.txt`** tutte le mutazioni stabilizzanti ($\Delta\Delta G < 0\text{ kcal/mol}$).

### Quali posizioni scansioniamo adesso?
Lo script `scripts/04_ankyrin_consensus_and_targets.py` (Passo 3.1) genera automaticamente il file di configurazione **`03_mutations/position_scan_hotspots.cfg`** che scansiona in tutti i 20 amminoacidi (`d`) i **principali hotspot strutturali e di consenso** già individuati nei 4 Ankyrin Repeats (tutti rigorosamente **fuori** dalla Functional Mask di 24 residui!):
1. **`LA73d` (`LEU 73` in `ANK3`, pos. 6 del repeat):** Dove tutti gli Ankyrin Repeats canonici hanno l'Istidina del motivo firma **`TPLH`**, `2MYO` ha una Leucina (`T70-P71-L72-L73`). Vediamo cosa dice FoldX se ripristiniamo **`L73H`** (o altri residui)!
2. **`KA94d` (`LYS 94` in `ANK3`, pos. 27 del repeat):** Nel motivo conservato `GADVNA` (dove la posizione 27 è una tasca idrofobica occupata da **Valina `V`** in `ANK1` e **Isoleucina `I`** in `ANK2`), `ANK3` ha infilato una **Lisina carica `LYS 94`** (`GADKTV`) con Ramachandran anomalo! Vediamo quanto stabilizza mutarla in **`K94V`**, **`K94I`** o **`K94L`**!
3. **`PA64d` (`PRO 64` in `ANK2`, pos. 30 del repeat):** Ha l'angolo di Ramachandran più distorto di tutta la proteina ($\psi = +121.45^\circ, \phi = -27.74^\circ$ nel `REMARK 500`).
4. **`LA105d` (`LEU 105` in `ANK4`, C-cap):** L'unico residuo che distingue il Ratto (`LEU 105`) dall'Uomo (`PHE 105` in `3AAA`)! Testiamo se la variante umana **`L105F`** (o un altro idrofobico di core) è più stabile!
5. **Più tutti gli altri siti chiave di inter-repeat/consenso e Ramachandran non mascherati** (`EA26d`, `RA30d`, `EA33d`, `GA35d`, `KA37d`, `LA48d`, `AA59d`, `HA68d`, `VA81d`, `SA82d`, `AA92d`, `DA93d`, `TA95d`, `KA97d`, `GA101d`, `TA103d`, `TA108d`, `NA110d`).
*(E domani mattina, se la MD a 400 K farà emergere qualche altro residuo con $\langle Q_i \rangle < 0.6$ non ancora incluso, basterà aggiungerlo in 10 secondi!)*

### 💻 Comando 3.3 — Lanciare `PositionScan` sui candidati (~5–8 minuti)
```bash
foldx -f 03_mutations/position_scan_hotspots.cfg > 03_mutations/position_scan.log
```
Una volta terminato, lancia il nostro script Python di analisi dello screening per classificare tutte le mutazioni e generare il grafico vettoriale `.pdf`:
```bash
python3 scripts/05_parse_foldx_screening.py
```
- **Output numerico:** `03_mutations/foldx_top_mutations.txt` (classifica ordinata dal $\Delta\Delta G$ più negativo/stabilizzante al meno stabilizzante).
- **Output grafico vettoriale:** `06_presentation/fig7_foldx_position_scan_heatmap.pdf`.

---

## 🌐 PASSO 3.4 — Lanciare i Predittori Ortogonali sul Cloud (0% CPU locale, lavorano mentre dormi!)

Mentre il tuo PC fa girare NAMD stanotte, possiamo far lavorare gratis i server accademici del **Pilastro 3** (`Lecture05`) caricando il file **`03_mutations/2MYO_Repair.pdb`**:

1. **ThermoMPNN (Graph Neural Network per $\Delta\Delta G$, *Dieckhaus & Kuhlman, PNAS 2024* — `Lecture05` Slide 65–69):**
   - Apri il Colab ufficiale: [ThermoMPNN Colab](https://colab.research.google.com/github/Kuhlman-Lab/ThermoMPNN/blob/main/ThermoMPNN.ipynb)
   - Carica `03_mutations/2MYO_Repair.pdb` (Chain `A`).
   - In **15 secondi** su GPU gratuita di Google, ThermoMPNN calcola il $\Delta\Delta G$ predetto per **tutte le $118 \times 20 = 2.360$ mutazioni possibili** di `2MYO` e ti fa scaricare il file `.csv`! Salvalo in `03_mutations/thermompnn_2MYO.csv`.
2. **FireProt 2.0 (Rosetta + FoldX + Back-to-Consensus Evolutivo — *Musil et al., NAR 2024*):**
   - Vai su [https://loschmidt.chemi.muni.cz/fireprot/](https://loschmidt.chemi.muni.cz/fireprot/), carica `03_mutations/2MYO_Repair.pdb` (o inserisci il PDB ID `2MYO`, Chain `A`) e avvia il job automatico.
3. **DynaMut2 (Normal Mode Analysis & Graph Signatures):**
   - Vai su [https://biosig.lab.uq.edu.au/dynamut2/](https://biosig.lab.uq.edu.au/dynamut2/prediction_mm) e testa la lista delle nostre migliori mutazioni emerse da FoldX + Consenso Ankyrin.

---

## 🏗️ PASSO 3.5 (Da eseguire domani mattina dopo la MD) — Matrice di Consenso Finale & `BuildModel`

Domani mattina, appena la MD del Wild-Type avrà prodotto `02_wt_baseline/unstable_residues_report.txt`:
1. Incroceremo la colonna **$\langle Q_i \rangle_{400\text{K}}$** (dalla MD) con **$\Delta\Delta G_{\text{FoldX}}$** (Passo 3.3), **Ankyrin Consensus** (Passo 3.1) e **ThermoMPNN/FireProt** (Passo 3.4).
2. Sceglieremo ufficialmente:
   - **`Mut1`** (Prima mutazione singola)
   - **`Mut2`** (Seconda mutazione singola)
   - **`Double`** (`Mut1 + Mut2`, il Doppio Mutante combinato)
3. Scriveremo le 3 righe in **`03_mutations/individual_list.txt`** (con il punto e virgola `;` finale obbligatorio, come in `Tutorials/foldx_test/build_model/individual_list.txt`, es.:
   ```text
   LA73H;
   KA94V;
   LA73H,KA94V;
   ```
   ) e lanceremo **`foldx -f 03_mutations/build_model.cfg`** (`numberOfRuns=3`) per ottenere i 3 PDB mutati e la scomposizione energetica (`Dif_2MYO_Repair.fxout`)!

---

## 📸 MOMENTO FOTO PER LE SLIDE (Fase 3 — Da salvare in `06_presentation/`)

Dalla Fase 3 otterremo **3 figure chiave per le Slide 5 e 6**:
1. **`06_presentation/fig6_ankyrin_consensus_alignment.pdf`** — Allineamento vettoriale dei 4 Ankyrin Repeats di `2MYO` contro il consenso universale di Mosavi et al. (2002), con evidenziati i residui mascherati (`3AAA`) e i residui "frustrati" (`LEU 73`, `LYS 94`, `PRO 64`, `LEU 105`).
2. **`06_presentation/fig7_foldx_position_scan_heatmap.pdf`** — Grafico vettoriale delle migliori mutazioni identificate da FoldX `PositionScan` ($\Delta\Delta G < -1.0\text{ kcal/mol}$) con scomposizione energetica.
3. **Foto 3D VMD (`scripts/viz_04_mutations_zoom.tcl`)** — Zoom affiancato in VMD (*Before vs After*) della tasca locale nel Wild-Type (`2MYO_Repair.pdb`) contro il Mutante generato da `BuildModel`, mostrando i nuovi contatti idrofobici o legami a idrogeno formati dalla mutazione!

---

## ✅ Checklist Operativa della Fase 3 (`ROAD3`)

- [ ] **Passo 3.1 (Adesso):** Eseguire `python3 scripts/04_ankyrin_consensus_and_targets.py` per generare il report del consenso Ankyrin, il grafico `.pdf` (`fig6`) e il file di configurazione `03_mutations/position_scan_hotspots.cfg`.
- [ ] **Passo 3.2 (Adesso):** Eseguire `foldx -f 03_mutations/repair.cfg` per riparare `2MYO.pdb` $\rightarrow$ `03_mutations/2MYO_Repair.pdb`.
- [ ] **Passo 3.3 (Adesso):** Eseguire `foldx -f 03_mutations/position_scan_hotspots.cfg` e poi `python3 scripts/05_parse_foldx_screening.py` per ottenere la classifica $\Delta\Delta G_{\text{FoldX}}$ e il grafico `.pdf` (`fig7`).
- [ ] **Passo 3.4 (Adesso / Stasera):** Sottomettere `03_mutations/2MYO_Repair.pdb` a ThermoMPNN / FireProt 2.0.
- [ ] **Passo 3.5 (Domani dopo la MD):** Incrociare con $\langle Q_i \rangle_{400\text{K}}$ della MD, scegliere `Mut1`, `Mut2` e `Double`, ed eseguire `foldx -f 03_mutations/build_model.cfg`.
