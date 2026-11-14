# 🌐 01_GUIDA_OPERATIVA_ALESSIA.md — Guida Completa allo Screening Mutazionale Online (Cloud / AI / Rosetta / NMA)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Proteina Target:** Myotrophin (*Rattus norvegicus*, PDB ID: `2MYO`, 118 amminoacidi, Catena `A`)  
**Presentazione Finale:** Mercoledì 16 Dicembre (max 10 minuti)

---

## 🎯 1. Cosa Stiamo Cercando e Perché il Tuo Lavoro è Fondamentale

### L'obiettivo del progetto in 3 righe
Vogliamo ingegnerizzare la **Myotrophin (`2MYO`)** per renderla **più termostabile** (cioè capace di resistere ad alte temperature senza perdere la propria struttura 3D, $\Delta\Delta G < 0$) individuando:
- **2 mutazioni singole indipendenti (`Mut1` e `Mut2`)**
- **1 doppio mutante combinato (`Double = Mut1 + Mut2`)**

### Perché non ci basta FoldX da solo (e perché entra in gioco la tua parte)?
Oleg sta facendo girare sul suo PC:
1. La **Dinamica Molecolare (MD) a `300 K` e `400 K`** per vedere quali amminoacidi del Wild-Type si "sbottonano" davvero quando l'acqua bolle;
2. Lo **screening energetico locale con FoldX 5** su tutte le posizioni modificabili della proteina.

Tuttavia, **FoldX da solo ha un limite noto** (`Lecture05` del Prof. Gautieri): tiene lo scheletro (*backbone*) della proteina quasi fermo e usa una funzione empirica che a volte produce "falsi positivi" (ad esempio premia amminoacidi troppo grossi o rigidi che poi in acqua calda rompono il fold).  
Per fare un vero **lavoro di ricerca** (e avere una parte solida e autonoma da raccontare nelle tue slide il 16 Dicembre!), affianchiamo a FoldX **3 predittori ortogonali sul Cloud** che usi tu dal browser (**0 installazioni sul tuo PC!**):

| Tool Online (Tuo Compito) | Su cosa si basa (`Lecture05`) | Cosa cerca nella nostra proteina `2MYO` |
| :--- | :--- | :--- |
| **1. FireProt 2.0** *(Server Accademico)* | **Evoluzione (MSA) + Rosetta + FoldX** | Confronta `2MYO` con centinaia di proteine omologhe naturali e ricalcola l'energia con **Rosetta** (molto più severo di FoldX sugli ingombri sterici nel core idrofobico). |
| **2. ThermoMPNN-D** *(Google Colab GPU)* | **Deep Learning (Graph Neural Network)** | Rete neurale addestrata su centinaia di migliaia di temperature di fusione sperimentali (*Megascale dataset*): in 15 secondi valuta tutte le $118 \times 20 = 2.360$ mutazioni singole e le coppie di doppi mutanti. |
| **3. DynaMut2** *(Server Accademico)* | **Normal Mode Analysis (Vibrazioni NMA)** | Controlla che le mutazioni candidate non introducano tensioni meccaniche o perdite di flessibilità anomale nelle vibrazioni della proteina. |
| **4. PROSS** *(Opzionale Bonus)* | **Filogenesi PSSM + Rosetta Design** | Secondo server indipendente (Weizmann Institute) basato su Rosetta per avere un'ulteriore conferma incrociata. |

---

## 🛑 2. La Regola d'Oro Biologica: Cosa Possiamo Toccare e Cosa è VIETATO Toccare

La Myotrophin (`2MYO`) è fatta da **4 mattoncini ripetuti** chiamati **Ankyrin Repeats** (`ANK1`, `ANK2`, `ANK3`, `ANK4`, 118 amminoacidi in tutto):
- **`ANK1`:** residui `1–34` (N-cap)
- **`ANK2`:** residui `35–67` (repeat centrale)
- **`ANK3`:** residui `68–100` (repeat centrale)
- **`ANK4`:** residui `101–118` (C-cap)

Nella cellula, la faccia concava della Myotrophin si lega all'**Actin-Capping Protein (CP)** per regolare il citoscheletro del muscolo cardiaco (struttura cristallografica del complesso: PDB `3AAA`).  
Se mutassimo un amminoacido che tocca l'Actin-Capping Protein, la proteina diventerebbe magari più stabile al calore, ma **smetterebbe di funzionare biologicamente** (`Laboratory_Project_2026.pdf`, Slide 7 & 19)!

### ❌ I 24 Residui VIETATI ("Functional Mask" calcolata a $< 4.5\text{ Å}$ su `3AAA`)
Qualsiasi mutazione proposta dai siti web su questi **24 numeri di residuo** va **SCARTATA SEMPRE**:
> **`4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36, 41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100`**

*(Escludiamo inoltre le prime e ultime 3 posizioni delle code libere **`1, 2, 3`** e **`116, 117, 118`**, che fluttuano nell'acqua).*

### ✅ Cosa Cerchiamo Davvero (Le 88 Posizioni Buone dello Scaffold!)
Cerchiamo mutazioni **fortemente stabilizzanti** sulle restanti **88 posizioni interne**, con particolare attenzione a:
1. **Ripristino del Consenso delle Ankyrin Repeats** (*Mosavi et al., PNAS 2002*): punti dove `2MYO` del ratto ha un amminoacido "anomalo" rispetto al mattoncino Ankyrin perfetto (es. **`A59G`**, **`H68G`/`H68N`**, **`T103L`/`T103F`**, **`L105F`** che è la variante naturale umana in `3AAA`, **`K37L`/`K37I`**, **`E26T`**, **`L73H`**, **`K94V`/`K94I`**).
2. **Rilascio di Tensioni Geometriche nei Loop Inter-Repeat**: punti con angoli di Ramachandran distorti nell'NMR (es. **`P64V`**, **`D60L`/`D60A`**, **`D93G`**, **`N29G`**).
3. **Rinforzo del Core Idrofobico del C-cap (`ANK4`) e di `ANK3`**: dove la proteina tende ad aprirsi ad alta temperatura (es. **`T108L`**, **`V81L`**, **`T103F`**).

---

## 📂 3. Quali File Usare (Già Pronti nella Cartella `03_mutations/` di Git!)

Dopo aver fatto `git pull`, dentro la cartella **`03_mutations/`** trovi i **2 soli file di input** che ti servono per tutti i siti web:

1. **`03_mutations/2MYO_Repair.pdb`**
   - **Cos'è:** È la struttura 3D di `2MYO` (Catena **`A`**, residui `1–118`) **già riparata con `FoldX RepairPDB`**.
   - **⚠️ IMPORTANTE:** Carica **SEMPRE `2MYO_Repair.pdb`** (e **NON** `2MYO.pdb` grezzo né il codice `2MYO` diretto se il sito ti permette di caricare il tuo file PDB). *Motivo scientifico:* `2MYO.pdb` viene da misure NMR e contiene piccoli scontri tra atomi che farebbero sballare i calcoli di Rosetta e DynaMut2; `2MYO_Repair.pdb` ha le catene laterali già rilassate ed è identico a quello usato da Oleg in FoldX!
2. **`03_mutations/dynamut2_mutations_list.txt`**
   - **Cos'è:** È il file di testo già formattato con l'elenco delle migliori mutazioni singole emerse da FoldX + Consenso Ankyrin (una per riga, formato `A D60L`, `A T108L`, ecc.). Ti serve **solo** per il Sito 3 (**DynaMut2**), perché DynaMut2 non fa la scansione cieca di tutta la proteina ma richiede una lista di mutazioni da testare.

---

## 🚀 4. Guida Passo-Passo ai 3 Siti Web (Ordine Consigliato)

> **💡 Consiglio sui tempi:** Lancia **per primo FireProt 2.0** (impieghi 3 minuti a compilarlo, ma poi il loro server ci mette da 1 a qualche ora per fare i calcoli con Rosetta e ti manda una mail). Subito dopo fai **ThermoMPNN-D** (impiega 15 secondi!) e **DynaMut2** (impiega 5–10 minuti).

---

### 🔥 SITO 1 — FireProt 2.0 (Evoluzione MSA + Rosetta + FoldX)
* **Link diretto:** [https://loschmidt.chemi.muni.cz/fireprot/](https://loschmidt.chemi.muni.cz/fireprot/)
* **Cosa fa:** Allinea la nostra sequenza contro centinaia di proteine simili in natura (*Back-to-Consensus*) e testa le mutazioni sia con **FoldX** sia con **Rosetta** (che valuta benissimo l'impaccamento idrofobico del core).
* **Passaggi esatti:**
  1. Apri il link e clicca su **Run FireProt** / **Submit a new job**.
  2. Nella scelta della struttura, seleziona **Upload PDB file** e carica il file **`03_mutations/2MYO_Repair.pdb`**.
  3. Quando ti mostra la catena proteica, seleziona la **Chain `A`** (118 residui).
  4. *(Opzionale ma utilissimo se il form te lo mostra nelle impostazioni avanzate "Fixed / Excluded residues"):* puoi indicare come **posizioni fisse da NON mutare** i 24 residui della nostra Functional Mask (`4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36, 41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100`). Se lasci le impostazioni standard non preoccuparti: il nostro script Python scarterà comunque in automatico i 24 residui mascherati!
  5. Inserisci la tua **email istituzionale (`@mail.polimi.it`)** o personale per ricevere il link appena il calcolo Rosetta finisce, e clicca **Submit**.
* **Cosa scaricare quando finisce:**
  - Nella pagina dei risultati di FireProt vedrai le mutazioni divise in *Energy-based* (FoldX/Rosetta) e *Evolution-based* (Consenso).
  - Clicca sul pulsante di **Download / Export CSV** (o scarica la tabella completa delle mutazioni valutate) e salva il file dentro la cartella `03_mutations/` chiamandolo:
    👉 **`03_mutations/fireprot_results.csv`**
  - Fai anche **1 o 2 screenshot** della vista 3D/tabella di FireProt (saranno perfetti per la tua slide!).

---

### 🤖 SITO 2 — ThermoMPNN-D su Google Colab (Intelligenza Artificiale — 15 secondi!)
* **Link diretto al Notebook Colab ufficiale (Kuhlman Lab):**  
  [https://colab.research.google.com/github/Kuhlman-Lab/ThermoMPNN-D/blob/main/ThermoMPNN-D.ipynb](https://colab.research.google.com/github/Kuhlman-Lab/ThermoMPNN-D/blob/main/ThermoMPNN-D.ipynb)
* **Cosa fa:** È la Graph Neural Network mostrata a lezione (`Lecture05`, Slide 65–69). Nella versione aggiornata **ThermoMPNN-D**, calcola sia tutte le **$118 \times 20 = 2.360$ mutazioni singole** (`single` mode), sia le **coppie di doppi mutanti** (`epistatic` / `additive` mode)!
* **Passaggi esatti su Google Colab:**
  1. Apri il link con il tuo account Google.
  2. In alto a destra (o dal menu `Runtime` $\rightarrow$ `Change runtime type`), verifica che sia selezionata la GPU gratuita (**`T4 GPU`**).
  3. Esegui la **Cella 1 (Setup / Installazione)** cliccando sul tasto `▶` a sinistra della cella (ci mette ~1 minuto a installare PyTorch e i pesi della rete).
  4. Nella **Cella di Upload / Input PDB**, carica dal tuo PC il file **`03_mutations/2MYO_Repair.pdb`** (Chain `A`).
  5. **Run 1 — Tutte le Mutazioni Singole (`mode = single`):**
     - Imposta la modalità su **`single`** (e se c'è il parametro `threshold`, mettilo a `100` o lascialo a `-0.5` per vedere tutte quelle stabilizzanti).
     - Premi `▶`: in **15 secondi** calcola i $\Delta\Delta G$ e genera il file `.csv`!
     - Scarica il file `.csv` delle mutazioni singole e salvalo in `03_mutations/` con il nome esatto:
       👉 **`03_mutations/thermompnn_2MYO.csv`**
  6. **Run 2 (Bonus Ricerca) — Mutazioni Doppie Epistatiche (`mode = epistatic`):**
     - Se nel notebook c'è la cella o l'opzione per **`epistatic`** (o `double`), lanciala sullo stesso `2MYO_Repair.pdb` e scarica anche quel CSV salvandolo come:
       👉 **`03_mutations/thermompnn_double_2MYO.csv`**
     - *Perché è super figo per la presentazione:* ci dice direttamente quali coppie di mutazioni lavorano in sinergia secondo l'AI!

---

### 🌊 SITO 3 — DynaMut2 (Analisi Vibrazionale & Normal Mode Analysis)
* **Link diretto:** [https://biosig.lab.uq.edu.au/dynamut2/](https://biosig.lab.uq.edu.au/dynamut2/)
* **Cosa fa:** A differenza di FoldX e ThermoMPNN, **DynaMut2** modella la proteina come una rete di molle armoniche (*Normal Mode Analysis*) e calcola se la mutazione rende la struttura più stabile senza irrigidire o destabilizzare i modi normali di vibrazione.
* **⚠️ ATTENZIONE ALLA PAGINA GIUSTA DEL MENU:**
  1. Vai su [https://biosig.lab.uq.edu.au/dynamut2/](https://biosig.lab.uq.edu.au/dynamut2/).
  2. Nel menu in alto clicca su **`Run` $\rightarrow$ `Single Mutation`** (oppure clicca sul box **Single Mutation** nella home page).
     *(Nota: **NON** andare su `Multiple Mutations` per questo file, perché `dynamut2_mutations_list.txt` contiene una mutazione singola per riga!).*
  3. Dentro la pagina **Single Mutation**, cerca il riquadro **Mutation List** (modalità batch per testare tante mutazioni insieme):
     - **Wild-type structure (Upload PDB file):** carica **`03_mutations/2MYO_Repair.pdb`**.
     - **Mutation list (Upload file):** carica **`03_mutations/dynamut2_mutations_list.txt`**.
     - Inserisci la tua email (opzionale ma comodo) e clicca **Submit**.
* **Cosa scaricare quando finisce (~5–10 minuti):**
  - Nella pagina dei risultati clicca su **Download CSV** (o tabella risultati) e salva il file in `03_mutations/` con il nome esatto:
    👉 **`03_mutations/dynamut2_results.csv`**

---

### 🧪 SITO 4 (Opzionale Bonus se hai 3 minuti extra) — PROSS (Weizmann Institute)
* **Link diretto:** [https://pross.weizmann.ac.il/](https://pross.weizmann.ac.il/)
* **Cosa fa:** È uno dei server più famosi al mondo per la termostabilità proteica (usa allineamento omologhi + filtro energetico **Rosetta**).
* **Nota:** Richiede obbligatoriamente l'iscrizione/invio con **email accademica (`@mail.polimi.it`)**. Se vuoi lanciarlo, carica `03_mutations/2MYO_Repair.pdb`, Chain `A`, e salva i risultati che ti arrivano via mail in `03_mutations/pross_results.txt`.

---

## ⚠️ 5. Come Leggere i Numeri che Escono (ATTENZIONE AI SEGNI!)

Quando guardi le tabelle che hai scaricato per farti un'idea di quali mutazioni vincono, fai **molta attenzione alla convenzione del segno di $\Delta\Delta G$**, perché **DynaMut2 usa il segno opposto rispetto agli altri!**:

| Programma | Quando una mutazione è **STABILIZZANTE** (Buona ✅) | Quando è **DESTABILIZZANTE** (Cattiva ❌) |
| :--- | :--- | :--- |
| **FoldX 5** *(Oleg)* | **$\Delta\Delta G < 0$** *(es. `-2.53 kcal/mol` = ottima!)* | $\Delta\Delta G > 0$ |
| **ThermoMPNN-D** *(Tu)* | **$\Delta\Delta G < 0$** *(es. `-1.20 kcal/mol` = ottima!)* | $\Delta\Delta G > 0$ |
| **FireProt 2.0 / Rosetta** *(Tu)* | **$\Delta\Delta G < 0$** *(es. `-1.80 kcal/mol` = ottima!)* | $\Delta\Delta G > 0$ |
| **DynaMut2** *(Tu — ⚠️ SEGNO OPPOSTO!)* | **$\Delta\Delta G > 0$ (POSITIVO!)** *(es. `+0.85 kcal/mol` = stabilizzante!)* | $\Delta\Delta G < 0$ (Negativo = destabilizzante!) |

*(Non preoccuparti per l'unione finale: lo script Python [`scripts/06_build_consensus_matrix.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/06_build_consensus_matrix.py) inverte già da solo il segno di `dynamut2_results.csv` così nella tabella finale tutti i numeri negativi significano "stabilizzante"!).*

---

## ✅ 6. Checklist Finale per Alessia (Cosa Caricare su Git)

Appena hai scaricato i file CSV dai siti, mettili nella cartella **`03_mutations/`** con questi nomi precisi:
- [ ] `03_mutations/thermompnn_2MYO.csv` *(da ThermoMPNN-D Colab — singole)*
- [ ] `03_mutations/thermompnn_double_2MYO.csv` *(opzionale, da ThermoMPNN-D Colab — doppie)*
- [ ] `03_mutations/dynamut2_results.csv` *(da DynaMut2)*
- [ ] `03_mutations/fireprot_results.csv` *(da FireProt 2.0 appena arriva la mail)*

Poi carica tutto sul repository Git con:
```bash
git add 03_mutations/*.csv
git commit -m "Add cloud mutation screening results (ThermoMPNN, DynaMut2, FireProt)"
git push
```
*(In alternativa, se preferisci non fare push da terminale, manda direttamente i 3 file `.csv` a Oleg su WhatsApp/Telegram e li mette lui in `03_mutations/`!).*
