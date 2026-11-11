# 💡 discovery.md — Registro delle Scoperte Scientifiche, Strutturali e HPC

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Target:** Myotrophin (*Rattus norvegicus*, PDB ID: `2MYO`)  

Questo file raccoglie tutte le **"chicche" scientifiche, biofisiche e computazionali (HPC)** scoperte e verificate sperimentalmente durante lo svolgimento del progetto. Ogni punto qui documentato è materiale prezioso da valorizzare nelle **Slide della Presentazione Finale (16 Dicembre)** per distinguere il nostro lavoro.

---

## 🧬 1. Scoperte Strutturali e Chimico-Fisiche (Fase 1)

### 1.1 La Rete di Legami a Idrogeno Inter-Repeat di `HIS 40` (`HSD`) e `HIS 67` (`HSE`)
- **Il problema standard:** Nei tutorial base di Dinamica Molecolare (es. `Tutorials/charmm36/psfgen.tcl`), tutte le istidine (`HIS`) vengono rinominate alla cieca con un'unica istruzione globale (`pdbalias residue HIS HSE` oppure `HSD`).
- **Cosa abbiamo scoperto misurando le distanze in `2MYO.pdb` (`scripts/01_inspect_2myo.tcl`):**
  1. **`HIS 40` (Istidina conservata del motivo firma `TPLH` nel Repeat 2):**
     - L'azoto imidazolico **`ND1`** si trova a soli **`2.44 Å`** dall'ossigeno idrossilico **`OG` di `SER 74`** (situata nel Repeat 3 adiacente).
     - Assegnando **`HSD`** (protone `HD1` legato a `ND1`), `HIS 40` dona un fortissimo legame a idrogeno inter-repeat che tiene cuciti insieme `ANK2` e `ANK3`. Se avessimo usato `HSE`, questo ponte idrogeno chiave sarebbe andato perso!
  2. **`HIS 67` (Situata nel $\beta$-hairpin loop tra Repeat 2 e Repeat 3):**
     - Nel modello NMR originale `2MYO.pdb`, X-PLOR aveva posizionato di default il protone su `ND1` (`HD1`) su tutte le istidine.
     - Tuttavia, la nostra ispezione geometrica entro $3.5\text{ Å}$ ha rivelato che:
       - A **`3.24 Å` da `ND1`** si trova il gruppo ammonico carico positivamente **`NZ` ($-\text{NH}_3^+$) di `LYS 37`**, che può **solo donare** legami a idrogeno (avere un protone anche su `ND1` genererebbe una repulsione elettrostatica/sterica diretta `+` contro `+`!).
       - A **`2.40 Å` da `NE2`** si trova l'ossigeno carbonilico del backbone **`O` (`C=O`) di `GLY 35`**, che è un **puro accettore** di legami a idrogeno!
     - Assegnando **`HSE`** (protone `HE2` su `NE2` e doppietto libero su `ND1`), `HIS 67` forma un **doppio ponte a idrogeno cooperativo**: accetta su `ND1` da `LYS 37` e dona da `NE2-HE2` a `GLY 35`, cucendo il loop tra `ANK1`, `ANK2` e `ANK3`!
- **Valore per le Slide (Slide 3):** Dimostra un'attenzione chimico-fisica di livello da ricerca (*Zero Ipotesi*), evitando di introdurre instabilità artificiali nel Wild-Type prima ancora di iniziare la MD.

---

### 1.2 Confronto Evoluzionistico Ratto (`2MYO`) vs Uomo (`3AAA`): 99.2% di Identità (Differisce solo `L105F`!)
- Confrontando riga per riga il record `SEQRES` della Myotrophin di Ratto (`2MYO`, `Chain A`, 118 aa, NMR in soluzione) con la Myotrophin Umana nel complesso cristallografico con l'Actin-Capping Protein (`3AAA`, `Chain C`, Raggi X a $2.20\text{ Å}$):
  - **117 amminoacidi su 118 sono identici ($99.2\%$ di identità di sequenza)!**
  - L'**unica sostituzione naturale** tra Ratto e Uomo su tutta la proteina si trova nel C-cap (`ANK4`), al residuo **`105`**: **`LEU 105`** nel Ratto (`2MYO`) contro **`PHE 105`** nell'Uomo (`3AAA`).
  - L'**RMSD globale sui $C_\alpha$ (`resid 2–118`)** tra la struttura NMR libera in soluzione (`2MYO`) e la struttura cristallografica legata al suo target biologico (`3AAA`) è di appena **`2.365 Å`**, dimostrando che la Myotrophin agisce come uno scaffold rigido pre-organizzato!

---

### 1.3 Mappatura Oggettiva della "Functional Mask" (24 Residui Intoccabili)
- Sovrapponendo `2MYO` a `3AAA` e selezionando tutti i residui di Myotrophin con atomi pesanti a distanza **$< 4.5\text{ Å}$** dalle subunità $\alpha$ e $\beta$ dell'Actin-Capping Protein (`Chain A` e `Chain B`), abbiamo identificato esattamente **24 residui funzionali "read-only"**:
  - **`ANK1`:** `LYS 4, GLU 5, MET 7, TRP 8, ALA 9, LYS 11, ASN 12, ASP 14, ASP 16, GLU 17, GLY 34, ARG 36`
  - **`ANK2`:** `TYR 41, ASP 44, CYS 45, GLN 47, LYS 66, HIS 67, ILE 69`
  - **`ANK3`:** `TYR 77, GLU 78, HIS 80, PRO 99, ASP 100`
- Questa maschera include al 100% la triade calda nota in letteratura (**`TRP 8, ASP 44, GLU 78`**, *Takeda et al., PLoS Biology 2010*) e lascia esattamente **94 residui di "telaio" (scaffold)** su cui progettare le mutazioni termostabilizzanti senza compromettere la funzione biologica.

---

### 1.4 "Frustrazioni Esecutive" di `2MYO` rispetto al Consenso Universale degli Ankyrin Repeats
- Allineando i 4 moduli (`ANK1–ANK4`) di `2MYO` contro il consenso termostabile universale degli Ankyrin Repeats di 33 aa (`G-X-TPLHLAA-X-X-G-H-L-E-I-V-E-V-L-L-K-X-G-A-D-V-N-A-X-D`, *Mosavi et al., PNAS 2002*), abbiamo scoperto 3 anomalie strutturali naturali nella Myotrophin:
  1. **Nel Repeat 3 (`ANK3`, pos. 73):** Il motivo firma `TPLH` è mutato in **`T70 - P71 - L72 - L73`** — manca l'Istidina conservata in posizione 73, sostituita da una Leucina apolare (`LEU 73`), privando l'interfaccia `ANK3–ANK4` del classico ponte idrogeno inter-repeat!
  2. **Nel Repeat 3 (`ANK3`, pos. 94–95):** Il motivo conservato `GADVNA` (presente quasi intatto in `ANK1` come `GEDVNR` e in `ANK2` come `GADINA`) in `ANK3` diventa **`G91 - A92 - D93 - K94 - T95 - V96`**: in posizione **`94`**, dove il consenso richiede una Valina/Isoleucina idrofobica (`V`/`I`), c'è una **Lisina carica positivamente (`LYS 94`)** (segnalata anche tra i torsionali anomali nel `REMARK 500` di `2MYO.pdb`!).
  3. **Nel Repeat 2 (`ANK2`, pos. 64):** Nel loop `GADINAPD` è presente **`PRO 64`**, che nel `REMARK 500` di `2MYO.pdb` presenta angoli di Ramachandran fortemente anomali ($\psi = +121.45^\circ, \phi = -27.74^\circ$).

---

## ⚡ 2. Scoperte Computazionali e di Performance HPC (Fase 2)

### 2.1 Contesa delle Unità Vettoriali AVX-512 in SMT: `+p8` (Core Fisici) batte `+p16` (SMT) del **+40%**!
- **Hardware di test:** CPU AMD 8 Core Fisici / 16 Thread Logici (2-way SMT) con supporto **AVX-512** (`AVX512F`, `AVX512CD`) e GPU integrata AMD Radeon 860M, 24 GB RAM.
- **Risultato sperimentale sul sistema `2myo_wb_ionized` (`21.461` atomi, PME, $\Delta t = 2\text{ fs}$):**
  - Con **`namd3 +p16`** (usando tutti i 16 thread logici in HyperThreading/SMT):
    - Tempo medio per step: **`0.0184 - 0.0191 s/step`**
    - Throughput: **`9.01 - 9.38 ns/day`** (`0.111 days/ns`)
  - Con **`namd3 +p8`** (usando esclusivamente gli **8 core fisici**):
    - Tempo medio per step: **`0.0133 - 0.0136 s/step`**
    - Throughput: **`12.51 - 12.68 ns/day`** (`0.077 days/ns`)
- **Spiegazione Architetturale HPC (da portare nella Slide 3 del Computational Budget!):**
  - Nei kernel compute-bound e fortemente vettorizzati di NAMD 3 (calcolo delle forze non-bonded di Lennard-Jones in **AVX-512** + FFT 3D per il **Particle Mesh Ewald** su una scatola compatta di `21.461` atomi divisa in 36 patch), lanciare 2 thread SMT sullo stesso core fisico provoca **contesa sulle pipeline vettoriali FMA/AVX-512, aumento del thrashing della cache L1/L2 e raddoppio dell'overhead di sincronizzazione/migrazione del Load Balancer Charm++** (`800 compute objects` su 16 PE invece di 8 PE).
  - Fissare il numero di worker agli **8 core fisici (`+p8`)** regala uno **speedup netto del $+40\%$** (`da 9.0 a 12.7 ns/day`) e riduce il tempo della campagna notturna di ben 3 ore!
