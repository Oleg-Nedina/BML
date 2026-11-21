# 🎯 PUNTO 3 DELLA SCALETTA DEL PROFESSORE — "3. Mutation hypotheses and screening: Identify positions and types of mutations, and shortlist them (e.g., using FoldX)"

**Riferimento Ufficiale nella Scaletta del Corso:**  
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 3** (*Learning Objective 2: Mutation design — Explain how substitutions could improve local interactions*), **Slide 7** (*Protein Engineering: Given the nearly infinite potential sequence space $20^{100} \approx 1.26 \times 10^{130}$, a computational filter is essential*), **Slide 12** (*Local contacts can guide candidate selection*), **Slide 15** (*Project workflow — Step 3: Mutation hypotheses and screening. Identify positions and types of mutations, and shortlist them*), **Slide 16–19** (*FoldX: Energetic screening $\Delta\Delta G = \Delta G_{\text{mut}} - \Delta G_{\text{WT}}$; FoldX workflow: `RepairPDB` $\rightarrow$ `PositionScan` $\rightarrow$ `BuildModel` using matched WT references $\rightarrow$ `Selection for MD` combining energetic prediction with structural plausibility and functional constraints*).

---

## 1. 📖 Introduzione al Problema e Riferimenti alle Slide di Teoria del Corso

Per una proteina di $118$ amminoacidi come la Myotrophin (`2MYO`), lo spazio delle sequenze possibili è pari a $20^{118} \approx 3.3 \times 10^{153}$ combinazioni (`Laboratory_Project_2026.pdf`, Slide 7). Anche limitandosi alle sole mutazioni puntiformi singole sulle `88` posizioni interne dello scaffold non-mascherato, esistono **$88 \times 19 = 1.672$ mutanti singoli** e oltre **$1.38\text{ milioni}$ di doppi mutanti**!  
Testarli tutti con la Dinamica Molecolare o con calcoli alchemici di integrazione termodinamica (`Lecture05 - Beyond Conventional MD.pdf`, Slide 35–40: *Alchemical transformations & Mutation free energy calculations*) è computazionalmente impossibile. Occorre quindi un **filtro di screening energetico e strutturale rapido** (`Laboratory_Project_2026.pdf`, Slide 7 & 16–19; `Lecture05`, Slide 12–16, 48–68):

1. **Termodinamica della Stabilità Proteica ($\Delta G_{\text{fold}}$ e $\Delta\Delta G$, `Lecture05` Slide 12–16; `Laboratory_Project_2026.pdf` Slide 5 & 17):**
   - La stabilità termodinamica di una proteina è governata dall'energia libera di Gibbs di folding $\Delta G = \Delta H - T\Delta S$ (`Lecture05`, Slide 12–14).
   - L'effetto di una mutazione puntiforme viene quantificato dalla variazione di energia libera di folding (`Laboratory_Project_2026.pdf`, Slide 17):
     $$\Delta\Delta G = \Delta G_{\text{mut}} - \Delta G_{\text{WT}} \quad (\text{kcal/mol})$$
     Nella convenzione di **FoldX**, **Rosetta/FireProt** e **ThermoMPNN**:
     - **$\Delta\Delta G < 0\text{ kcal/mol}$ (Negativo):** mutazione **stabilizzante** (sposta l'equilibrio verso lo stato folded e innalza la $T_m$, specialmente per $\Delta\Delta G < -1.0\text{ kcal/mol}$);
     - **$\Delta\Delta G > 0\text{ kcal/mol}$ (Positivo):** mutazione **destabilizzante**;
     - **$\Delta\Delta G \approx 0\text{ kcal/mol}$:** effetto neutro o entro il rumore del predittore ($\pm 0.5\text{ kcal/mol}$).
2. **Il Workflow Ufficiale di FoldX (`Laboratory_Project_2026.pdf`, Slide 16–19; `Tutorials/foldx_test/`):**
   - FoldX utilizza una funzione di energia empirica lineare calibrata sperimentalmente:
     $$\Delta G_{\text{FoldX}} = W_{\text{vdw}}\Delta G_{\text{vdw}} + W_{\text{solvH}}\Delta G_{\text{solvH}} + W_{\text{solvP}}\Delta G_{\text{solvP}} + \Delta G_{\text{hbond}} + \Delta G_{\text{el}} + \Delta G_{\text{clash}} + T\Delta S_{\text{mc}} + T\Delta S_{\text{sc}}$$
   - Il protocollo rigoroso prescritto nella **Slide 19 di `Laboratory_Project_2026.pdf`** si articola in 4 fasi:
     1. **`RepairPDB`:** ottimizza i rotameri delle catene laterali del PDB nativo per eliminare gli scontri sterici di Van der Waals prima dello screening.
     2. **`PositionScan`:** muta sistematicamente ciascun residuo selezionato in tutti i 20 amminoacidi (`d`), rilassando i residui vicini per identificare le sostituzioni con $\Delta\Delta G < 0$.
     3. **`BuildModel`:** costruisce i modelli atomici 3D dei mutanti singoli e combinati finalisti eseguendo più run indipendenti (`numberOfRuns=5`) e generando per ciascuna run un **riferimento Wild-Type accoppiato (*matched WT reference*, `WT_2MYO_Repair_*.pdb`)** per calcolare $\Delta\Delta G_{\text{BuildModel}}$ senza bias conformazionali.
     4. **`Selection for MD`:** combina la predizione energetica con la **plausibilità strutturale** e i **vincoli funzionali** (*"FoldX predicts an energetic effect. MD tests structural behaviour under chosen conditions"*, Slide 17 & 19).
3. **L'Upgrade "Research-Grade" dalla `Lecture05` (Slide 48–68) e dal Paper ScALDO (2025): Consenso Ortogonale Multi-Tool:**
   - Come spiegato nella `Lecture05` (Slide 48–59 su **Rosetta**: *Rotamers, Side-chain packing & Monte Carlo sampling*; Slide 60–68 su **Inverse Folding & ProteinMPNN**: *Backbone geometry as a graph*) e nel paper del Journal Club del corso (*2025 - Substrate Promiscuity Engineering of ScALDO*), un singolo predittore a backbone semi-rigido come FoldX può generare **falsi positivi** (ad es. premiando l'inserimento di Proline in mezzo a un'$\alpha$-elica o la rimozione di un accettore H-bond non standard).
   - Per eliminare ogni falso positivo, affianchiamo a FoldX e alla MD a `400 K` tre pilastri ortogonali:
     - **Consenso Strutturale degli Ankyrin Repeats** (*Mosavi et al., PNAS 2002*; *Lowe & Itzhaki, PNAS 2007*);
     - **FireProt 2.0** (*Musil et al., NAR 2024* — combina *Back-to-Consensus* filogenetico MSA `isBTC` con filtri energetici FoldX/Rosetta, `Lecture05` Slide 48–59);
     - **ThermoMPNN-D** (*Dieckhaus & Kuhlman, PNAS 2024* — Graph Neural Network basata su **ProteinMPNN**, `Lecture05` Slide 66–68, addestrata sul *Megascale dataset* di $\Delta\Delta G$ sperimentali per mutazioni singole e coppie epistatiche doppie).

---

## 2. 🎯 Cosa Cercavamo (Gli Obiettivi del Punto 3)

1. **Riparazione Energetica del Modello NMR (`RepairPDB`):** Eliminare i clash sterici e torsionali dell'NMR grezzo `2MYO.pdb` per ottenere una struttura di riferimento `2MYO_Repair.pdb` energeticamente rilassata.
2. **Mappatura delle "Frustrazioni Esecutive" rispetto al Consenso Ankyrin:** Individuare dove i 4 Ankyrin Repeats naturali della Myotrophin deviano dalla sequenza consenso termostabile di 33 residui (*Mosavi et al., 2002*) o presentano angoli di Ramachandran $(\phi, \psi)$ fuori regione (`REMARK 500` di `2MYO.pdb`).
3. **Copertura Energetica Totale (`88 su 88` posizioni dello scaffold) con `FoldX PositionScan`:** Calcolare il $\Delta\Delta G_{\text{FoldX}}$ e la scomposizione fisico-chimica esatta ($\Delta E_{\text{VdW+SolvH}}$, $\Delta E_{\text{Clash}}$, $\Delta E_{\text{HBond}}$, $\Delta E_{\text{Elec}}$, $\Delta E_{\text{Entropy}}$) per **tutte le $88 \times 20 = 1.760$ varianti possibili** dello scaffold interno non-mascherato.
4. **Convergenza a 5 Pilastri (*Orthogonal Consensus Funnel*):** Unire i dati di **MD WT a `400 K` ($\langle Q_i \rangle_{400\text{K}}$)** + **FoldX 5 (`88/88`)** + **Consenso Ankyrin** + **FireProt 2.0 (`fireprot_results.csv`)** + **ThermoMPNN-D (`thermompnn_2MYO.csv` & `thermompnn_double_2MYO.csv`)** per selezionare i migliori mutanti singoli complementari (`Mut1`, `Mut2`) e il corrispondente Doppio Mutante (`Double = Mut1 + Mut2`), spiegando chimicamente *come* ciascuna sostituzione migliora le interazioni locali (*Learning Objective 2*, Slide 3).

---

## 3. 🛠️ Come L'Abbiamo Cercato (Metodologia, Script e Dati Unificati)

1. **Passo 3.1 — Allineamento Consenso Ankyrin & Outlier di Ramachandran ([`scripts/04_ankyrin_consensus_and_targets.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/04_ankyrin_consensus_and_targets.py)):**
   - Allineati i 4 repeat (`ANK1` 2–34, `ANK2` 35–67, `ANK3` 68–100, `ANK4` 101–118) alle 33 posizioni canoniche di *Mosavi et al. (2002)*, incrociandoli con i 35 outlier di Ramachandran del `REMARK 500` e mascherando i 24 residui funzionali di `3AAA`.
   - Generati [`03_mutations/ankyrin_consensus_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ankyrin_consensus_report.txt) e il grafico vettoriale [`06_presentation/fig6_ankyrin_consensus_alignment.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig6_ankyrin_consensus_alignment.pdf).
2. **Passo 3.2 — `FoldX RepairPDB` ([`03_mutations/repair.cfg`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/repair.cfg)):**
   - Eseguito `foldx -f 03_mutations/repair.cfg` su `2MYO.pdb` generando [`03_mutations/2MYO_Repair.pdb`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/2MYO_Repair.pdb) e il trace di ottimizzazione [`03_mutations/2MYO_Repair.fxout`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/2MYO_Repair.fxout).
3. **Passo 3.3 — `FoldX PositionScan` Parallelo a 8 Core sul 100% dello Scaffold (`88/88` posizioni, [`scripts/04b_run_remaining_foldx_parallel.sh`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/04b_run_remaining_foldx_parallel.sh) + [`scripts/05_parse_foldx_screening.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/05_parse_foldx_screening.py)):**
   - Tolti i `24` residui della Functional Mask e i `6` residui delle code terminali libere (`1–3` e `116–118`), abbiamo scansionato con `PositionScan` (`d` = 20 amminoacidi) **tutte le 88 posizioni interne modificabili della proteina** (`43` hotspot nel primo batch + `45` posizioni rimanenti distribuite su 8 worker paralleli), generando `88` file `03_mutations/energies_<res>_2MYO_Repair.txt`, la classifica [`03_mutations/foldx_top_mutations.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/foldx_top_mutations.txt) e il grafico vettoriale [`06_presentation/fig7_foldx_position_scan_heatmap.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig7_foldx_position_scan_heatmap.pdf).
4. **Passo 3.4 — Screening Cloud Evolutivo & AI (Lavoro di Alessia su `2MYO_Repair.pdb`):**
   - **FireProt 2.0:** sottomesso `2MYO_Repair.pdb` (`Chain A`) ed esportato [`03_mutations/fireprot_results.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/fireprot_results.csv) (`2.242` mutazioni con `foldx`, `conservation_score`, `correlation_score`, `majority`, `ratio`, `bfactor`, `isBTC`, `isHighEnergy`, `isHighCombined`).
   - **ThermoMPNN-D (Colab GPU):** eseguito sia in modalità `single` ([`03_mutations/thermompnn_2MYO.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/thermompnn_2MYO.csv), `2.242` mutazioni) sia in modalità `double/epistatic` ([`03_mutations/thermompnn_double_2MYO.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/thermompnn_double_2MYO.csv), `46.440` coppie).
5. **Passo 3.5 — Fusione nella Matrice di Consenso a 5 Pilastri ([`scripts/06_build_consensus_matrix.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/06_build_consensus_matrix.py)):**
   - Applicazione automatica delle **6 Regole** di [`MUTATION_SELECTION.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/MUTATION_SELECTION.md), calcolo del numero di pilastri concordi (`Pillars / 5`), generazione di [`03_mutations/consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt), di [`06_presentation/fig8_consensus_energy_breakdown.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8_consensus_energy_breakdown.pdf) e dell'analisi monografica [`03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md).

---

## 4. 💡 Cosa Abbiamo Trovato e Considerazioni Critiche (Risultati Quantitativi)

### 4.1 Effetto di `FoldX RepairPDB` sulla Struttura NMR ([`2MYO_Repair.fxout`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/2MYO_Repair.fxout))
- **Energia Totale FoldX prima di `RepairPDB` (`2MYO.pdb`):** **`+319.97 kcal/mol`** (`VdW Clash = 69.13 kcal/mol`, `Torsional Clash = 86.95 kcal/mol`, `H-Bond = -33.12 kcal/mol`).
- **Energia Totale FoldX dopo `RepairPDB` (`2MYO_Repair.pdb`):** **`+170.53 kcal/mol`** (**$\Delta E_{\text{repair}} = -149.44\text{ kcal/mol}$**!).
  - I clash di Van der Waals sono crollati dell'**80.4%** (da `69.13` a **`13.53 kcal/mol`**), i clash torsionali del **60.1%** (da `86.95` a **`34.67 kcal/mol`**) e l'energia dei legami a idrogeno è migliorata da `-33.12` a **`-56.81 kcal/mol`** (`-35.73` backbone + `-21.08` sidechain).

### 4.2 I Vincitori della Matrice di Consenso a 5 Pilastri ([`consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt))
Dall'incrocio oggettivo dei 5 pilastri emergono **3 mutazioni singole d'eccellenza situate nei 3 moduli che mostrano vulnerabilità a 400 K (`ANK2`, `ANK3`, `ANK4`)**:

| Rank / Candidato | Repeat & Posizione Strutturale | Pilastri Concordi | $\Delta\Delta G_{\text{FoldX}}$ (kcal/mol) | Scomposizione FoldX ($\Delta E_{\text{VdW+SolvH}}$ / $\Delta E_{\text{Clash}}$) | ThermoMPNN-D (kcal/mol) | FireProt 2.0 ($\Delta\Delta G$ & Flags) | MD WT $\langle Q_i \rangle_{400\text{K}}$ & RMSF | Come Migliora le Interazioni Locali (*Learning Objective 2*, Slide 3) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **🥇 `K37L`**<br>*(Mutante Primario `ANK2`)* | **`ANK2`** (pos. 3 del repeat, ingresso motivo `TPLH` `37–40`) | **5 / 5** | **`-2.876`** | `-0.255` / `-0.127`<br>*(+ solvatazione polare favorevole `-1.84`)* | **`-0.359`**<br>*(Sito #1 dell'AI!)* | **`-0.89`** | **`0.4845`** (`UNSTABLE`)<br>`RMSF = 2.96 Å` | **1)** Rimuove la carica positiva di `LYS 37` adiacente a `ARG 36` (eliminando la repulsione elettrostatica `+`/`+`).<br>**2)** Abolisce il ponte salino parassita `LYS37-NZ -> ASP60-OD1/OD2` (`47.6%` occupancy a 300 K) che deformava il loop `GADVNA` di `ANK2`.<br>**3)** La catena isobutilica della Leucina sigilla idrofobicamente il cap N-terminale dell'elica interna di `ANK2` (il modulo più instabile a 400 K, $Q_{\text{ANK2}} = 0.586$). |
| **🥈 `T103L`**<br>*(Mutante Primario `ANK4` C-Cap)* | **`ANK4`** (pos. 3/5 del C-cap, ingresso elica `102–110`) | **5 / 5** | **`-1.670`**<br>*(o `T103F`: `-2.103`)* | **`-3.476`** / `+1.053` | **`-0.116`**<br>*(o `T103F`: `-0.313`)* | **`-1.18`** | **`0.5075`** (`UNSTABLE`)<br>`RMSF = 3.16 Å`<br>*(accanto a `L102: 0.4942`)* | **1)** Ripristina la **Leucina (`L`) del consenso Ankyrin** (*Mosavi et al., 2002*) all'ingresso dell'elica del C-cap.<br>**2)** Guadagna ben **`-3.48 kcal/mol` di energia idrofobica e di Van der Waals**, ancorando l'inizio dell'elica C-cap (`LEU 102 – THR 103`, che a 400 K perdeva il legame H `ALA107 -> LEU102`) contro il core di `ANK3`.<br>**3)** Stabilizza direttamente il **nucleo di folding primario** C-terminale della Myotrophin (*Lowe & Itzhaki, PNAS 2007*). |
| **🥉 `H68G`**<br>*(Mutante #1 Evolutivo FireProt in `ANK3`)* | **`ANK3`** (pos. 1 del repeat, apice $\beta$-turn `ANK2–ANK3`) | **4 / 5**<br>*(Score #1: `-6.98`)* | **`-1.736`** | `+0.440` / **`-0.670`**<br>*(+ $\Delta S_{\text{mc}} = -1.29$)* | `+0.029`<br>*(neutra)* | **`-2.33`**<br>**(`BTC, MAJ, RAT, HighC`)** | **`0.7139`** (`WEAK`)<br>`RMSF = 3.63 Å`<br>*(accanto a `H67: 0.4336`)* | **1)** È la **#1 mutazione *Back-to-Consensus* di FireProt 2.0**.<br>**2)** Ripristina la **Glicina (`G`) universale in posizione 1** del modulo Ankyrin (`G-X-TPLHLAA`, *Mosavi et al., 2002*).<br>**3)** Elimina l'outlier di Ramachandran nel quadrante sinistrorso ($\phi = +38.28^\circ, \psi = +25.85^\circ$, `REMARK 500` riga 146), accessibile senza strain sterico solo alla Glicina! |

*(Inoltre, nella Tabella Master abbiamo documentato anche i candidati di supporto **`D93G`** [`4/5 pilastri`, `FP = -2.27 HighE/HighC`, `FoldX = -1.564`, $\langle Q_{93} \rangle = 0.5427$], **`T108A`/`T108L`** in `ANK4`, e la coppia locale sinergica **`A59G + D60L`** in `ANK2` con epistasi ThermoMPNN-D $\Delta\Delta G_{\text{double}} = -0.756\text{ kcal/mol}$ e $\sum\Delta\Delta G_{\text{FoldX}} = -4.945\text{ kcal/mol}$!).*

### 4.3 Considerazioni Critiche: Perché la Nostra Selezione Batte uno Screening Automatico Cieco
Come documentato nel dettaglio in [`03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md), il nostro imbuto a 5 pilastri ha intercettato e **scartato 3 classi di falsi positivi** che avrebbero rovinato la proteina in MD:
1. **La "trappola" `S74A` (rottura del ponte H inter-repeat `HIS 40 -> SER 74`):** Sia FireProt (`isBTC = true`, `-0.55`), sia FoldX (`-2.395`), sia ThermoMPNN (`-0.17`) proponevano `S74A`. Solo grazie alla nostra ispezione del Punto 1 sapevamo che `SER 74 : OG` accetta l'H-bond inter-repeat ($2.44\text{ Å}$) da `HIS 40 (HSD)`!
2. **Gli *helix-breakers* `E52P` (`FP = -2.28`) e `T103G` (`FP = -2.60`):** Proposti da FireProt/FoldX ma situati dentro $\alpha$-eliche (`HELIX 4` e C-cap helix) e smascherati sia dalla nostra Regola 3 sia da ThermoMPNN-D (`+0.44` e `+0.69 kcal/mol`).
3. **Le 10 mutazioni nella Functional Mask (`D16W` e `I69G/M/R/N/S/H/K/L/D`):** Proposte da FireProt ma situate a $< 4.5\text{ Å}$ dall'Actin-Capping Protein in `3AAA`.

---

## 5. ✅ Risultati di `FoldX BuildModel` (`n=5` Run con *Matched WT References*, Slide 19) e Chiusura al 100% del Punto 3

Eseguendo [`scripts/07_run_foldx_buildmodel.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/07_run_foldx_buildmodel.py) (`numberOfRuns=5` contro i riferimenti Wild-Type accoppiati `WT_2MYO_Repair_*.pdb`, report completo in [`03_mutations/buildmodel_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/buildmodel_epistasis_report.txt) e grafico vettoriale [`06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf)), abbiamo ottenuto la validazione termodinamica definitiva pre-MD (*Livello 1 / Gate 1 Superato*):

| Variante (`BuildModel`, $n=5$) | Tipo | $\Delta\Delta G_{\text{BuildModel}}$ ($\mu \pm \sigma$, kcal/mol) | $\Delta\Delta G_{\text{Best}}$ (kcal/mol) | Somma Singole ($\sum \Delta\Delta G_i$) | Accoppiamento Epistatico ($\Delta\Delta\Delta G_{\text{coupling}}$) | Verdetto Termodinamico & Strutturale (*Gate 1*) |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **`K37L`** (`ANK2`, 5/5★) | Singola | **`-2.438 ± 0.870`** | **`-3.324`** | — | — | ✅ **Mutante Singolo 1 (`Mut1`):** forte guadagno di solvatazione polare (`-2.02`) e idrofobica (`-0.68`) ed entropia (`-1.07`). |
| **`T103L`** (`ANK4`, 5/5★) | Singola | **`-1.150 ± 0.143`** | **`-1.250`** | — | — | ✅ **Mutante Singolo 2 (`Mut2`):** convergenza strettissima ($\sigma = 0.14$) e guadagno idrofobico enorme nel core C-cap (**$\Delta E_{\text{VdW+SolvH}} = -4.027\text{ kcal/mol}$**). |
| **`H68G`** (`ANK3`, #1 BTC) | Singola | **`-1.730 ± 0.009`** | **`-1.747`** | — | — | ✅ **Mutante Singolo 3 (`Mut3 / BTC`):** convergenza deterministica ($\sigma = 0.009$), rilascio di clash di backbone (`-0.60`) ed entropia di mainchain (`-1.35`). |
| **`K37L + T103L`** (`ANK2 + ANK4`) | **Doppia A** | **`-3.652 ± 0.656`** | **`-4.766`** | `-3.588` | **`-0.064`** | 🏆 **DOPPIO MUTANTE PRIMARIO (`Double`):** **100% Additivo** con leggera sinergia (`-0.064 kcal/mol`) e zero interferenza sterica tra i due poli del solenoide! |
| **`H68G + T103L`** (`ANK3 + ANK4`) | Doppia C | **`-2.888 ± 0.418`** | **`-3.592`** | `-2.880` | **`-0.008`** | ✅ **100% Additivo** ($\Delta\Delta\Delta G = -0.008\text{ kcal/mol}$). |
| **`K37L + H68G`** (`ANK2 + ANK3`) | Doppia B | `-2.889 ± 0.076` | `-2.982` | `-4.168` | `+1.278` | ⚠️ **Sub-Additivo:** poiché `K37` e `H68` sono vicini sulla stessa cresta dei $\beta$-hairpin attorno a `HIS 67`, i loro gusci di rilassamento locale si sovrappongono parzialmente (`+1.28 kcal/mol`), dimostrando quanto sia stato prezioso eseguire `BuildModel` prima della MD! |
| **`A59G + D60L`** (`ANK2` locale) | Doppia E | **`-3.548 ± 0.676`** | **`-4.236`** | `-2.985` | **`-0.563`** | 🌟 **Sinergia Epistatica Positiva:** conferma in FoldX (`-0.563 kcal/mol` di sinergia) l'epistasi locale predetta dalla rete neurale ThermoMPNN-D (`A59G` sblocca la tensione sterica di `D60L`!). |
| **`K37L + H68G + T103L`** (`ANK2+3+4`) | **Tripla** | **`-5.032 ± 0.210`** | **`-5.392`** | `-5.318` | `+0.287` | 🚀 **Record Assoluto di Stabilità (`-5.03 kcal/mol`, $\sigma = 0.21$):** combina un "bullone" per ciascuno dei 3 repeat (`ANK2`, `ANK3`, `ANK4`). |

### Come Agiamo Adesso nel Punto 4 (Sistemi Confermati per le Simulazioni MD Comparative):
I conformeri 3D a più bassa energia sono già stati estratti in [`03_mutations/final_models_pdb/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/final_models_pdb):
1. **`WT`** (`2MYO` Wild-Type — Replica 1 a 400 K già completata)
2. **`Mut1 = K37L`** (`03_mutations/final_models_pdb/2myo_mut_K37L.pdb`, `ANK2`, 5/5 pilastri)
3. **`Mut2 = T103L`** (`03_mutations/final_models_pdb/2myo_mut_T103L.pdb`, `ANK4`, 5/5 pilastri)
4. **`Double = K37L + T103L`** (`03_mutations/final_models_pdb/2myo_mut_K37L_T103L.pdb`, `ANK2 + ANK4`, perfettamente additivo `-3.65 ± 0.66 kcal/mol`, `Best = -4.77 kcal/mol`)
*(più opzionalmente **`H68G`** / **`K37L + H68G + T103L`** già pronti in `final_models_pdb/`).*
