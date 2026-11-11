# 🗺️ ROAD1.md — Fase 1: Analisi Strutturale, Biofisica e Setup Topologico di Myotrophin (`2MYO`)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Target:** Myotrophin (*Rattus norvegicus*), PDB ID: [`2MYO`](https://www.rcsb.org/structure/2MYO)  
**Complesso Funzionale di Riferimento:** Actin Capping Protein (CP) + Myotrophin (V-1), PDB ID: [`3AAA`](https://www.rcsb.org/structure/3AAA)  
**Cartella di lavoro di questa fase:** `01_structures/`  
**Regole di Riferimento:** Tutte le operazioni seguono tassativamente quanto stabilito in [`AGENTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/AGENTS.md) (*Zero Ipotesi, Citazione Obbligatoria della Documentazione, Esecuzione da parte dell'Utente*).

---

## 📚 1. Fonti e Documentazione di Riferimento per la Fase 1

Ogni passaggio di questa guida è basato esclusivamente sulle seguenti fonti documentali e di letteratura:

1. **Linee Guida Ufficiali del Progetto:**  
   [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 8** (*Learning Objective 1: Protein Analysis — Fold, Secondary Structure, Relevant Regions*), **Slide 15** (*Workflow Step 1: Study the target protein*).
2. **Slide delle Lezioni del Corso (Prof. A. Gautieri):**
   - [`Lecture01 - Introduction.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture01%20-%20Introduction.pdf) — **Slide 24–38** (Struttura primaria/secondaria/terziaria, formato PDB, risoluzione e metodi sperimentali X-Ray vs NMR).
   - [`Lecture03 - Molecular Dynamics Setup.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture03%20-%20Molecular%20Dynamics%20Setup.pdf) — **Slide 18–25** (File `.psf` e `.pdb`, tool `psfgen`, rimozione idrogeni pre-esistenti, stati di protonazione delle istidine `HSD`/`HSE`/`HSP`, ponti disolfuro), **Slide 35–41** (Modelli di solvente esplicito `TIP3P` vs implicito `GBIS`, Periodic Boundary Conditions, padding minimo del box, neutralizzazione ed equilibrazione ionica con `autoionize`).
3. **Script e Force Field Ufficiali dei Tutorial del Corso:**
   - [`Tutorials/charmm36/psfgen.tcl`](file:///home/olly/UNI/BioMolEng2026/Tutorials/charmm36/psfgen.tcl) & [`Tutorials/nativecontacts/psfgen.tcl`](file:///home/olly/UNI/BioMolEng2026/Tutorials/nativecontacts/psfgen.tcl) (Procedura `psfgen` e selezione `"protein and noh"`).
   - [`Tutorials/charmm36/top_all36_prot.rtf`](file:///home/olly/UNI/BioMolEng2026/Tutorials/charmm36/top_all36_prot.rtf), [`par_all36m_prot.prm`](file:///home/olly/UNI/BioMolEng2026/Tutorials/charmm36/par_all36m_prot.prm) e [`toppar_water_ions.str`](file:///home/olly/UNI/BioMolEng2026/Tutorials/charmm36/toppar_water_ions.str) (Force field **CHARMM36m**, *Huang et al., Nature Methods 2017*).
   - [`VMDextensions.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/VMDextensions.tcl) (T. Giorgino, 2019 — funzioni `getFasta`, `prepareNativeContacts`, `helicity`).
4. **Dati Primari dei File PDB [`2MYO.pdb`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/2MYO.pdb) e [`3AAA_complex_ref.pdb`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/3AAA_complex_ref.pdb) & Letteratura Primaria:**
   - **Yang Y, Nanduri S, Sen S, Qin J (1998).** *The structural basis of ankyrin-like repeat function as revealed by the solution structure of myotrophin.* **Structure**, 6(5):619–626.
   - **Takeda S et al. (2010).** *Two distinct mechanisms for actin capping protein regulation—steric and allosteric inhibition.* **PLoS Biology**, 8(7):e1000416 (Struttura cristallografica `3AAA` a 2.20 Å del complesso Actin-Capping Protein + V-1/Myotrophin).
   - **Mosavi LK, Minor DL, Peng ZY (2002).** *Consensus-derived structural determinants of the ankyrin repeat motif.* **PNAS**, 99(25):16029–16034.
   - **Lowe AR, Itzhaki LS (2007).** *Rational redesign of the folding pathway of a modular protein.* **PNAS**, 104(8):2679–2684.

---

## 🎯 2. Obiettivi Operativi della Fase 1 (`ROAD1`)

Al termine di questo documento avrete completato e verificato 6 punti fondamentali:

1. **Passo 1.1 — Ispezione del file PDB grezzo (`2MYO.pdb`):** Verifica sperimentale dei record `HEADER`, `EXPDTA`, `REMARK 210` (condizioni NMR: $T = 297\text{ K}$, $\text{pH} = 6.2$), `REMARK 350` (stato oligomerico), `DBREF`, `SEQRES` e `HELIX`.
2. **Passo 1.2 — Verifica di Completezza, Cisteine e Stati di Protonazione delle Istidine in VMD:** Controllo oggettivo (zero ipotesi!) delle distanze $S_\gamma\text{--}S_\gamma$ delle 3 Cisteine (**`CYS 2`**, **`CYS 45`**, **`CYS 83`**) e della rete di legami a idrogeno attorno alle 4 Istidine (**`HIS 40`**, **`HIS 67`**, **`HIS 68`**, **`HIS 80`**) per scegliere rigorosamente tra `HSD`, `HSE` e `HSP`.
3. **Passo 1.3 — Mappatura dell'Architettura a 4 Ankyrin Repeats e Allineamento col Consenso Universale (`Mosavi et al., 2002`):** Divisione strutturale dei 4 moduli (`ANK1`–`ANK4`) e identificazione delle posizioni in cui `2MYO` devia dal consenso termostabile degli Ankyrin Repeats.
4. **Passo 1.4 — "Functional Masking" Oggettivo tramite il Complesso Cristallografico `3AAA`:** Sovrapposizione in VMD di `2MYO` (`Chain A`) su `3AAA` (`Chain C`) per estrarre la lista esatta dei residui di Myotrophin che contattano l'Actin-Capping Protein (`Chain A` e `Chain B` di `3AAA` entro $4.5\text{ Å}$). Questi residui costituiranno la **maschera funzionale intoccabile**.
5. **Passo 1.5 — Generazione della Topologia `CHARMM36m` (`psfgen`):** Rimozione degli idrogeni NMR (`protein and noh`), applicazione degli alias verificati al Passo 1.2 e creazione di `01_structures/2myo_wt.psf` e `01_structures/2myo_wt.pdb`.
6. **Passo 1.6 — Solvatazione Esplicita (`TIP3P`) e Ionizzazione Fisiologica (`0.15 M NaCl`):** Creazione del sistema finale solvatato e neutralizzato (`2myo_wb_ionized.psf` / `.pdb`) e misura esatta dei vettori della cella periodica (`cellBasisVector1/2/3` e `cellOrigin`) per NAMD.

---

## 🔬 PASSO 1.1 — Analisi Documentale del File `2MYO.pdb` (Cosa dice il file ufficiale)

Prima ancora di aprire la grafica 3D, leggiamo i metadati ufficiali contenuti nelle righe `1–200` di [`01_structures/raw/2MYO.pdb`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/2MYO.pdb#L1-L200) (che avete appena scaricato da RCSB PDB):

| Campo PDB (Righe in `2MYO.pdb`) | Contenuto Verificato nel File | Implicazione Operativa per il nostro Setup MD (`Lecture03`) |
| :--- | :--- | :--- |
| **`HEADER` / `TITLE`** (Righe 1–2) | `SOLUTION STRUCTURE OF MYOTROPHIN, NMR, MINIMIZED AVERAGE STRUCTURE` | È una struttura ottenuta in soluzione via **NMR** (media minimizzata rappresentativa di 90 conformeri calcolati, `REMARK 210` righe 66–67). Contiene **1 singolo modello** (non serve estrarre il `MODEL 1` da un ensemble). |
| **`COMPND` / `REMARK 350`** (Righe 3–5, 92–94) | `CHAIN: A`, `AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC` | L'unità biologica funzionale in soluzione è un **monomero a singola catena (`Chain A`)**. |
| **`REMARK 210`** (Righe 48–53) | `TEMPERATURE: 297 K`, `PH: 6.2`, `IONIC STRENGTH: 50 MM PHOSPHATE BUFFER`, `SAMPLE CONTENTS: H2O` | La struttura nativa è stata risolta in acqua a temperatura ambiente ($297\text{ K}$) e $\text{pH } 6.2$. |
| **`DBREF`** (Riga 163) | `DBREF 2MYO A 2 118 UNP P62775 MTPN_RAT 1 117` | **Attenzione critica alla numerazione:** Nel PDB `2MYO`, il residuo `1` è la Metionina iniziale (`MET 1`), mentre nella sequenza matura UniProt (`P62775`) la Metionina iniziale è scissa (`CYS 2` del PDB = residuo `1` in UniProt). **Nel progetto useremo sempre e solo la numerazione PDB (`1..118`)**, facendo attenzione quando confrontiamo con gli articoli che usano la numerazione senza Met1! |
| **`SEQRES`** (Righe 164–173) | `118` residui totali (`MET 1` $\rightarrow$ `GLN 118`) | La catena polipeptidica è completa da 1 a 118. |
| **`HELIX`** (Righe 174–181) | **8 segmenti elicoidali** annotati nel PDB:<ul><li>`HELIX 1`: `ASP 3` – `LYS 11` (Repeat 1, elica interna)</li><li>`HELIX 2`: `LEU 15` – `VAL 22` (Repeat 1, elica esterna)</li><li>`HELIX 3`: `HIS 40` – `CYS 45` (Repeat 2, elica interna)</li><li>`HELIX 4`: `ILE 50` – `LEU 56` (Repeat 2, elica esterna)</li><li>`HELIX 5`: `ASP 60` – `ALA 63` (Breve elica nel loop inter-repeat 2–3)</li><li>`HELIX 6`: `PRO 71` – `GLU 78` (Repeat 3, elica interna)</li><li>`HELIX 7`: `SER 82` – `LEU 88` (Repeat 3, elica esterna)</li><li>`HELIX 8`: `LYS 114` – `LEU 116` (Helix $3_{10}$ C-terminale; nota: nel paper di Yang et al. 1998 il Repeat 4 comprende anche l'elica `LEU 102` – `ASN 110`)</li></ul> | Ci fornisce la mappa ufficiale della struttura secondaria nativa con cui confronteremo la stabilità delle eliche durante la MD ad alta temperatura! |
| **`CRYST1`** (Riga 182) & **`REMARK 215`** (Righe 74–78) | `CRYST1 1.000 1.000 1.000 ... VALUES ON THESE RECORDS ARE MEANINGLESS` | Essendo una struttura NMR (non un cristallo), il record `CRYST1` è fittizio (`1.0 Å`). Le dimensioni del box periodico andranno ricalcolate dopo `solvate` con `measure minmax`. |
| **`ATOM`** (Righe 189–2004) | `1815` atomi totali (`ATOM 1` – `ATOM 1815`), inclusi **902 atomi di idrogeno NMR** (`H1`, `H2`, `H3`, `HA`, `HB2`, `HB3`...) e l'ossigeno C-terminale `OXT` (`ATOM 1807`). Zero molecole d'acqua (`HOH`) o eteroatomi (`HETATM`). | Come prescritto in `Lecture03` (Slide 21) e in `Tutorials/nativecontacts/psfgen.tcl`, gli idrogeni NMR (`H1`, `HB2`...) hanno nomi incompatibili con la topologia CHARMM e causano errori in `psfgen`. Vanno rimossi selezionando **`protein and noh`**. |

---

## 🖥️ PASSO 1.2 — Ispezione Chimico-Fisica in VMD: Completezza, Cisteine e Istidine (Zero Ipotesi!)

Ora verifichiamo direttamente dalle coordinate 3D di `2MYO.pdb`:
1. La sequenza FASTA continua e la carica netta della proteina.
2. La distanza tra gli atomi di zolfo (`SG`) delle **3 Cisteine (`CYS 2`, `CYS 45`, `CYS 83`)** per verificare se esistono ponti disolfuro (un ponte disolfuro richiede una distanza $d(S_\gamma\text{--}S_\gamma) \approx 2.05\text{ Å}$).
3. L'intorno di legami a idrogeno entro $3.5\text{ Å}$ degli azoti imidazolici **`ND1`** e **`NE2`** delle **4 Istidine (`HIS 40`, `HIS 67`, `HIS 68`, `HIS 80`)**, oltre alla posizione dell'idrogeno sperimentale NMR già presente nel file `2MYO.pdb`!
   - **Fatto documentale straordinario:** Poiché `2MYO.pdb` è una struttura NMR risolta a pH 6.2 con idrogeni espliciti, possiamo leggere direttamente nel file PDB quali atomi di idrogeno (`HD1` su `ND1` e/o `HE2` su `NE2`) sono presenti nelle 4 istidine e quali legami a idrogeno formano!

### 💻 Comando 1.2 (Da eseguire per Oleg su Linux o Alessia su Windows)

Posizionatevi nella cartella principale della repository `BML` ed eseguite questo comando di ispezione in VMD (potete lanciarlo direttamente dal terminale: apre VMD in modalità testuale, analizza `01_structures/raw/2MYO.pdb` e stampa il report completo):

#### 🐧 Su Linux (Bash) / 🪟 Su Windows (Git Bash):
```bash
vmd -dispdev text -eofexit << 'EOF'
mol new 01_structures/raw/2MYO.pdb

puts "\n=========================================="
puts "1. SEQUENZA FASTA E COMPLETEZZA CATENA"
puts "=========================================="
set prot [atomselect top "protein"]
puts "Numero totale atomi nel PDB: [$prot num]"
set heavy [atomselect top "protein and noh"]
puts "Numero atomi pesanti (noh):  [$heavy num]"
set ca [atomselect top "protein and name CA"]
puts "Numero residui (CA):         [$ca num]"
puts "Sequenza FASTA (da VMDextensions):"
puts [getFasta $prot]

puts "\n=========================================="
puts "2. VERIFICA DELLE 3 CISTEINE (DISTANZE SG-SG)"
puts "=========================================="
set cys_sg [atomselect top "resname CYS and name SG"]
set cys_res [$cys_sg get resid]
set cys_crd [$cys_sg get {x y z}]
for {set i 0} {$i < [llength $cys_res]} {incr i} {
    for {set j [expr {$i+1}]} {$j < [llength $cys_res]} {incr j} {
        set r1 [lindex $cys_res $i]
        set r2 [lindex $cys_res $j]
        set d [vecdist [lindex $cys_crd $i] [lindex $cys_crd $j]]
        puts [format "Distanza SG(CYS %d) -- SG(CYS %d) = %.2f A" $r1 $r2 $d]
    }
}

puts "\n=========================================="
puts "3. ISPEZIONE DELLE 4 ISTIDINE (IDROGENI NMR E CONTATTI < 3.5 A)"
puts "=========================================="
foreach hres {40 67 68 80} {
    set hsel [atomselect top "resid $hres"]
    set anames [$hsel get name]
    puts "\n--- HIS $hres ---"
    puts "Atomi presenti in 2MYO.pdb per HIS $hres: $anames"
    foreach ring_n {ND1 NE2} {
        set n_sel [atomselect top "resid $hres and name $ring_n"]
        set nbrs [atomselect top "(noh and not resid $hres) and within 3.5 of (resid $hres and name $ring_n)"]
        puts "  Vicini polari/pesanti entro 3.5 A da $ring_n:"
        foreach r [$nbrs get resname] id [$nbrs get resid] an [$nbrs get name] crd [$nbrs get {x y z}] {
            set dist [vecdist [lindex [$n_sel get {x y z}] 0] $crd]
            puts [format "    -> %s %d : %s (d = %.2f A)" $r $id $an $dist]
        }
    }
}
EOF
```

*(🪟 **Nota per Windows PowerShell:** Se Alessia usa PowerShell invece di Git Bash, può salvare i comandi tra `mol new ...` e l'ultimo `}` in un file oppure incollarli direttamente nella finestra **Extensions $\rightarrow$ Tk Console** di VMD!).*

> **👉 Fermatevi un istante dopo aver lanciato questo comando e guardate l'output insieme:**
> 1. Verificherete che non ci sono `/` nella sequenza FASTA (118 residui continui, 913 atomi pesanti + 902 idrogeni = 1815 atomi).
> 2. Vedrete le 3 distanze $S_\gamma\text{--}S_\gamma$ tra `CYS 2`, `CYS 45` e `CYS 83` (tutte $> 20\text{ Å}$, quindi **nessun ponte disolfuro**!).
> 3. Vedrete quali atomi di idrogeno sono legati a `ND1` (`HD1`) o `NE2` (`HE2`) nelle 4 istidine e quali legami a idrogeno formano con i residui vicini (in particolare vedrete il famoso legame a idrogeno strutturale di **`HIS 40`** con il loop dell'Ankyrin Repeat!).

---

## 🧬 PASSO 1.3 — Architettura dei 4 Ankyrin Repeats e Allineamento col Consenso (`Mosavi et al., 2002`)

Dalla letteratura primaria sulla struttura di `2MYO` (**Yang et al., Structure 1998**) e sul consenso universale delle proteine ad Ankyrin Repeat (**Mosavi, Minor & Peng, PNAS 2002, 99:16029–16034**), ogni modulo Ankyrin canonico è lungo **33 amminoacidi** (posizioni $1\dots 33$ del motivo ripetuto) e presenta:
- **Posizioni 1–3:** Connettore / $\beta$-turn (`G-X-T` nelle posizioni 2–4).
- **Posizioni 4–7 (Motivo firma `TPLH`):** La Treonina (`T4`) e l'Istidina (`H7`) formano una rete conservata di legami a idrogeno che stabilizza l'inizio della prima $\alpha$-elica e aggancia il $\beta$-hairpin loop del repeat adiacente.
- **Posizioni 5–12 ($\alpha$-elica 1, interna/concava):** `P-L-H-L-A-A-X-X`.
- **Posizioni 13–14 (Loop inter-elica):** `G-H` (la Glicina in pos. 13 consente l'inversione di direzione tra le due eliche antiparallele).
- **Posizioni 15–22 ($\alpha$-elica 2, esterna/convessa):** `L-E-I-V-E-V-L-L`.
- **Posizioni 25–33 (Loop lungo / $\beta$-hairpin inter-repeat):** `G-A-D-V-N-A-X-D-X`.

### Allineamento Esatto dei 4 Ankyrin Repeats di `2MYO` contro il Consenso Universale

Leggendo la sequenza ufficiale dal record `SEQRES` di `2MYO.pdb` (righe 164–173), ecco l'allineamento posizione per posizione dei 4 repeat di `2MYO` contro il consenso di Mosavi et al. (2002):

```text
Posizione nel Repeat (1..33):  01 02 03 [04 05 06 07 08 09 10 11] 12 13 14 [15 16 17 18 19 20 21 22] 23 24 25 26 27 28 29 30 31 32 33
Struttura Secondaria:          --beta-- [----- Elica 1 (Int) ----] --Loop-- [----- Elica 2 (Est) ----] -------- Loop / Beta-Turn -------
CONSENSO ANKYRIN (Mosavi 02):   D  G  x [ T  P  L  H  L  A  A  x]  x  G  H [ L  E  I  V  E  V  L  L]  K  x  G  A  D  V  N  A  x  D  x
----------------------------------------------------------------------------------------------------------------------------------------
ANK1 (N-cap, PDB res   1-34):   -  -  M [ C  D  K  E  F  M  W  A]  L  K  N [ G  D  L  D  E  V  K  D]  Y  V  A  K  G  E  D  V  N  R  T  L  E  G
ANK2 (Core,  PDB res  35-67):   G  R  K [ P  L  H  Y  A  A  D  C]  G  Q  L [ E  I  L  E  F  L  L  L]  K  G  A  D  I  N  A  P  D  K  H
ANK3 (Core,  PDB res  68-100):  H  I  T [ P  L  L  S  A  V  Y  E]  G  H  V [ S  C  V  K  L  L  L  S]  K  G  A  D  K  T  V  K  G  P  D
ANK4 (C-cap, PDB res 101-118):  G  L  T [ A  L  E  A  T  D  N  Q]  A  I  K [ A  L  L  Q  -  -  -  -]
```

> [!TIP]
> **Osservazioni Scientifiche Immediate dall'Allineamento (Fondamentali per la Fase 3 e per le Slide!):**
> 1. **In `ANK2` (res. 35–67):**
>    - Il motivo `TPLH` in `2MYO` è **`K37 - P38 - L39 - H40`** (manca la Treonina in pos. 37, sostituita da `LYS 37`).
>    - Nel loop `G-H` inter-elica (`pos. 13–14`), `2MYO` ha **`G46 - Q47`** (`GLN 47` al posto di `HIS`).
>    - Nel motivo `GADVNA` (`pos. 25–30`), `2MYO` ha **`G58 - A59 - D60 - I61 - N62 - A63 - P64 - D65`**: è quasi identico al consenso perfetto, tranne **`PRO 64`** (che infatti nel `REMARK 500` riga 145 di `2MYO.pdb` è segnalato come outlier di Ramachandran: $\psi = 121.45^\circ, \phi = -27.74^\circ$!).
> 2. **In `ANK3` (res. 68–100):**
>    - Il motivo `TPLH` in `2MYO` è **`T70 - P71 - L72 - L73`** (**manca l'Istidina conservata in pos. 73**, sostituita da una Leucina idrofobica `LEU 73`!). Questa è una deviazione famosissima studiata anche da *Lowe & Itzhaki (PNAS 2007)*!
>    - Nel motivo `GADVNA` (`pos. 25–30`, res. 91–96), `2MYO` ha **`G91 - A92 - D93 - K94 - T95 - V96`**: in posizione `94` (normalmente occupata da un residuo idrofobico ramificato `VAL`/`ILE` che si impacca nel core idrofobico tra repeat 3 e 4) c'è una **`LYS 94` carica positivamente**, e in posizione `95` (dove il consenso ha un'Asparagina `ASN` che forma H-bond col backbone inter-repeat) c'è **`THR 95`**!
> 3. **In `ANK4` (res. 101–118):**
>    - Ricordate il confronto che abbiamo fatto prima tra `2MYO` (Ratto) e `3AAA` (Uomo)? L'unica differenza su 118 amminoacidi tra la Myotrophin di ratto e quella umana è proprio in `ANK4` al residuo **`105`**: **`LEU 105`** nel ratto (`2MYO`) vs **`PHE 105`** nell'uomo (`3AAA`)!

---

## 🛡️ PASSO 1.4 — "Functional Masking": Identificazione Oggettiva dell'Interfaccia con Actin-Capping Protein (`3AAA`)

Come richiesto nella Slide 9 e Slide 15 di [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf), le mutazioni non devono alterare i siti funzionali della proteina.  
Invece di limitarci a citare i 3 residui principali (`TRP 8`, `ASP 44`, `GLU 78`), usiamo VMD per **sovrapporre strutturalmente `2MYO` (`Chain A`) sulla Myotrophin del complesso cristallografico `3AAA` (`Chain C`)** e misurare oggettivamente tutti i residui di Myotrophin che si trovano entro **$4.5\text{ Å}$** (soglia standard dei contatti atomici, *Bekker et al., 2019*) dalle subunità $\alpha$ (`Chain A`) e $\beta$ (`Chain B`) dell'Actin-Capping Protein!

### 💻 Comando 1.4 (Misura oggettiva dei residui d'interfaccia in `3AAA` e RMSD tra `2MYO` NMR e `3AAA` X-Ray)

Eseguite questo comando dal terminale nella cartella `BML`:

```bash
vmd -dispdev text -eofexit << 'EOF'
# Carica 3AAA (mol 0) e 2MYO (mol 1)
mol new 01_structures/raw/3AAA_complex_ref.pdb
mol new 01_structures/raw/2MYO.pdb

puts "\n=========================================="
puts "1. RMSD DEL BACKBONE TRA 2MYO (NMR) E 3AAA Chain C (X-RAY 2.20 A)"
puts "=========================================="
set sel_3aaa [atomselect 0 "chain C and resid 2 to 118 and name CA"]
set sel_2myo [atomselect 1 "chain A and resid 2 to 118 and name CA"]
set trans_mat [measure fit $sel_2myo $sel_3aaa]
$sel_2myo move $trans_mat
set rmsd_val [measure rmsd $sel_2myo $sel_3aaa]
puts [format "RMSD (CA resid 2-118) tra 2MYO (Ratto NMR) e 3AAA (Uomo X-Ray): %.3f A" $rmsd_val]

puts "\n=========================================="
puts "2. FUNCTIONAL MASK: RESIDUI DI MYOTROPHIN A CONTATTO (< 4.5 A) CON CAPPING PROTEIN"
puts "=========================================="
set interface_sel [atomselect 0 "chain C and noh and within 4.5 of (chain A B and noh)"]
set mask_resids [lsort -unique -integer [$interface_sel get resid]]
puts "Lista residui funzionali di Myotrophin a distanza < 4.5 A da Capping Protein:"
foreach r $mask_resids {
    set r_ca [atomselect 0 "chain C and resid $r and name CA"]
    puts [format "  -> Resid %3d (%s)" $r [$r_ca get resname]]
}
puts "\nStringa pronta per esclusione in VMD/Python (Functional Mask):"
puts "resid $mask_resids"
EOF
```

> **👉 Cosa otterrete da questo comando:**
> - Il valore esatto in Å dell'**RMSD tra la struttura NMR in soluzione (`2MYO`) e la struttura cristallografica legata a CP (`3AAA`)** (un dato eccellente da mettere nella **Slide 2** della presentazione!).
> - L'elenco completo e certificato di tutti i residui che toccano l'Actin-Capping Protein: questa lista costituisce la nostra **Functional Exclusion Mask** ufficiale. Nessuna mutazione verrà mai proposta su questi residui!

---

## ⚙️ PASSO 1.5 — Generazione della Topologia `CHARMM36m` con `psfgen`

Una volta ispezionate le 4 istidine al Passo 1.2, siamo pronti a costruire la topologia `.psf` e le coordinate `.pdb` del Wild-Type secondo la procedura vista in [`Lecture03`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture03%20-%20Molecular%20Dynamics%20Setup.pdf) (Slide 20–25) e in [`Tutorials/nativecontacts/psfgen.tcl`](file:///home/olly/UNI/BioMolEng2026/Tutorials/nativecontacts/psfgen.tcl):

1. **Rimozione degli idrogeni NMR pre-esistenti:** selezioniamo `"chain A and protein and noh"` e salviamo `01_structures/2myo_clean_noh.pdb`.
2. **Caricamento della topologia `CHARMM36m`:** usiamo il file locale condiviso nella repo `common/charmm36/top_all36_prot.rtf`.
3. **Alias per Isoleucina e Istidine:**
   - `pdbalias atom ILE CD1 CD` (converte il nome PDB standard del carbonio $\delta$ dell'isoleucina nel nome CHARMM `CD`, come in `Tutorials/charmm36/psfgen.tcl`).
   - Impostiamo lo stato di protonazione (`HSD` o `HSE`) per ciascuna delle 4 istidine (`40`, `67`, `68`, `80`) in base a quanto verificato nell'output del Comando 1.2!
4. **Ricostruzione delle coordinate degli idrogeni (`guesscoord`)** e salvataggio di `01_structures/2myo_wt.psf` e `01_structures/2myo_wt.pdb`.

---

## 💧 PASSO 1.6 — Solvatazione Esplicita (`TIP3P`) e Ionizzazione (`0.15 M NaCl`)

Seguendo [`Lecture03`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture03%20-%20Molecular%20Dynamics%20Setup.pdf) (Slide 36–39) e il paper di riferimento (*Bekker et al., 2019*):
- Usiamo il pacchetto `solvate` di VMD per immergere `2myo_wt` in una scatola d'acqua `TIP3P` con un padding di **$12.0\text{ Å}$** su ogni lato (`-t 12`), garantendo che la distanza minima tra immagini periodiche sia $\ge 24.0\text{ Å}$ (superiore al doppio del cutoff delle interazioni non-bonded di $12.0\text{ Å}$).
- Usiamo il pacchetto `autoionize` di VMD per neutralizzare la carica netta della Myotrophin e portare la concentrazione salina a **$0.15\text{ M NaCl}$** (`-sc 0.15`, forza ionica fisiologica).
- Infine, misuriamo con `measure minmax` e `measure center` le dimensioni esatte della scatola (`cellBasisVector1`, `cellBasisVector2`, `cellBasisVector3` e `cellOrigin`) da inserire nel file `.conf` di NAMD della Fase 2!

*(Risultati verificati dai Comandi 1.2 e 1.4: zero ponti disolfuro sulle 3 Cisteine `CYS 2, 45, 83`; tautomeri assegnati: `HIS 40 -> HSD`, `HIS 67 -> HSE`, `HIS 68 -> HSD`, `HIS 80 -> HSD`; RMSD `2MYO` vs `3AAA` = `2.365 Å`; Functional Mask di 24 residui: `resid 4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100`).*

### 💻 Comando 1.5 & 1.6 (Costruzione Topologia `CHARMM36m`, Solvatazione `TIP3P` e Ionizzazione `0.15 M NaCl`)

Dalla cartella principale `BML/`, eseguite lo script [`scripts/02_build_wt_charmm36m.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/02_build_wt_charmm36m.tcl):

```bash
vmd -dispdev text -e scripts/02_build_wt_charmm36m.tcl
```

---

## 📸 MOMENTO FOTO PER LE SLIDE (Fase 1 — Da salvare in `06_presentation/`)

Dalla Fase 1 potete già ricavare **3 immagini spettacolari** per le **Slide 2 e Slide 3** della presentazione finale!  
Ho preparato uno script automatico ([`scripts/viz_01_structure_and_mask.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_01_structure_and_mask.tcl)) che imposta in VMD lo sfondo bianco da articolo scientifico, toglie gli assi, applica l'ombreggiatura *Ambient Occlusion* e colora automaticamente tutti gli elementi:

### Come aprire la scena grafica pronta per gli screenshot:
```bash
vmd -e scripts/viz_01_structure_and_mask.tcl
```

Una volta aperta la finestra di VMD, andate su **`Graphics -> Representations`** (cliccando due volte su una riga nella lista potete accenderla/spegnerla in un attimo!) e salvate queste **3 figure** (usando **`File -> Render -> Snapshot`** oppure semplicemente uno screenshot ritagliato, salvandole in `06_presentation/`):

1. **Figura 1 (`06_presentation/fig1_ankyrin_repeats_and_mask.png`) — *Per la Slide 2 (Introduzione a `2MYO`)*:**
   - **Cosa lasciare acceso in `Representations`:** I 4 Ankyrin Repeats (`ANK1` Blu, `ANK2` Ciano, `ANK3` Arancione, `ANK4` Verde) + la Triade Funzionale (`TRP 8, ASP 44, GLU 78` in bastoncini rossi) + opzionalmente fate doppio click sulla riga `Surf` (Rappresentazione 5) per mostrare in rosso semitrasparente tutta la faccia concava dei 24 residui della *Functional Mask* che legano l'Actin-Capping Protein!
   - **Inquadratura:** Ruotate la proteina con il mouse (tasto sinistro) in modo che si vedano chiaramente le 4 "spire" del solenoide e i loop che sporgono da un lato.
2. **Figura 2 (`06_presentation/fig2_histidine_hbond_network.png`) — *Per la Slide 3 (Rigore nel Setup: Tautomeri `HIS 40` e `HIS 67`)*:**
   - **Cosa lasciare acceso:** I repeat centrali (`ANK1`, `ANK2`, `ANK3`) e l'ultima rappresentazione (`resid 35 37 40 67 74` in `Licorice`).
   - **Trucco VMD per disegnare i trattini delle "calamite" (H-bonds) e le distanze in Å sullo schermo:**
     - Premete il tasto **`2`** sulla tastiera (attiva lo strumento *Bond/Distance Measure* del mouse) e cliccate prima sull'azoto `ND1` di `HIS 40` e poi sull'ossigeno `OG` di `SER 74` $\rightarrow$ comparirà la linea tratteggiata con scritto **`2.44`**!
     - Poi cliccate su `NZ` di `LYS 37` e `ND1` di `HIS 67` (**`3.24`**) e su `NE2` di `HIS 67` e `O` di `GLY 35` (**`2.40`**)! *(Per tornare a ruotare normalmente col mouse basta premere il tasto **`r`** sulla tastiera).*
3. **Figura 3 (`06_presentation/fig3_solvated_box.png`) — *Per la Slide 3 (Dominio di Simulazione PBC)*:**
   - Dopo aver lanciato `scripts/02_build_wt_charmm36m.tcl`, potete aprire la scatola d'acqua con:
     ```bash
     vmd 01_structures/2myo_wb_ionized.psf 01_structures/2myo_wb_ionized.pdb
     ```
   - Impostando `protein` in `NewCartoon`, `ions` (`resname SOD CLA`) in `VDW` (sfere) e `water` in `Lines` (o `QuickSurf` trasparente) e scrivendo nella Tk Console `pbc box` per disegnare il cubo del dominio periodico!

---

## ✅ Checklist di Controllo della Fase 1 (`ROAD1`)

Prima di passare a `ROAD2.md` (Benchmark NAMD e Simulazione WT 300 K vs 400 K), verificheremo insieme di avere:
- [x] Eseguito il **Comando 1.2** e annotato le distanze delle 3 Cisteine (`CYS 2, 45, 83`) e i legami a idrogeno/protonazione delle 4 Istidine (`HIS 40 -> HSD, HIS 67 -> HSE, HIS 68 -> HSD, HIS 80 -> HSD`).
- [x] Eseguito il **Comando 1.4** e annotato l'RMSD tra `2MYO` e `3AAA` (`2.365 Å`) e la lista esatta dei 24 residui della **Functional Mask**.
- [ ] Generato `01_structures/2myo_wt.psf` e `01_structures/2myo_wt.pdb` e il sistema solvatato/ionizzato `01_structures/2myo_wb_ionized.psf` e `.pdb` lanciando `vmd -dispdev text -e scripts/02_build_wt_charmm36m.tcl`.
- [ ] Salvato le prime figure da VMD in `06_presentation/` con `vmd -e scripts/viz_01_structure_and_mask.tcl`.
