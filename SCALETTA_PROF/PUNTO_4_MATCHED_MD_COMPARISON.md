# ⚙️ PUNTO 4 DELLA SCALETTA DEL PROFESSORE — "4. Matched MD comparison: Simulate WT and variants with the same protocol and multiple velocity seeds"

**Riferimento Ufficiale nella Scaletta del Corso:**  
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 3** (*Learning Objective 3: Controlled computational comparison — Use consistent modelling and MD settings for WT and variants*), **Slide 10** (*Reference study: simulation protocol — multiple velocity seeds per system*), **Slide 15** (*Project workflow — Step 4: Matched MD comparison. Simulate WT and variants with the same protocol and multiple velocity seeds*), **Slide 20–21** (*MD simulations and analysis: Computational budget and pilot benchmark — Use the same protocol for WT and variants*).

---

## 1. 📖 Introduzione al Problema e Riferimenti alle Slide di Teoria del Corso

Perché non basta il $\Delta\Delta G < 0$ calcolato da FoldX, Rosetta o ThermoMPNN nel Punto 3 per dimostrare che una variante è più termostabile? E perché non basta una singola simulazione MD per variante?

1. **Predizione Energetica Statica vs Comportamento Dinamico in Acqua Calda (`Laboratory_Project_2026.pdf`, Slide 17 & 19; `Lecture04`, Slide 10–27; `Lecture05`, Slide 57):**
   - Come ribadito nella **Slide 17 di `Laboratory_Project_2026.pdf`** (*"FoldX predicts an energetic effect. MD tests structural behaviour under chosen conditions"*) e nella tabella comparativa della **`Lecture05` Slide 57** (*Molecular dynamics and Rosetta Monte Carlo*), i predittori di design operano a backbone quasi fisso e in solvente implicito.
   - Solo la Dinamica Molecolare all-atom in solvente esplicito (`TIP3P` + ioni a `400 K`) permette al backbone di rilassarsi, alle molecole d'acqua di competere per i legami a idrogeno intramolecolari e agli effetti cooperativi a lungo raggio tra gli Ankyrin Repeats di manifestarsi nel tempo (`Lecture04`, Slide 10–14).
2. **Il Principio del "Matched Comparison" (`Laboratory_Project_2026.pdf`, Slide 3, 15, 21; `Lecture04`, Slide 26–27):**
   - Nella **Slide 3 (*Controlled computational comparison*)** e nella **Slide 21 (*"Use the same protocol for WT and variants"*)**, il Prof. Gautieri stabilisce una regola metodologica aurea: affinché il confronto tra Wild-Type (`WT`) e mutanti (`Mut1`, `Mut2`, `Double`) sia scientificamente valido, **ogni singolo parametro del modello e del solutore deve essere mantenuto rigorosamente identico**:
     - Stesso force field (`CHARMM36m`), stesso modello d'acqua (`TIP3P`), stesso padding del box ($12.0\text{ Å}$), stessa forza ionica ($0.15\text{ M NaCl}$) e stessa neutralità elettrica ($q_{\text{tot}} = 0.0000\text{ e}$);
     - Stessi stati di protonazione dei residui non mutati (`HIS 40 -> HSD`, `HIS 67 -> HSE`, `HIS 68 -> HSD` [o `GLY` se mutata], `HIS 80 -> HSD`);
     - Stesso protocollo di minimizzazione (`1.000` step), stesso timestep ($\Delta t = 2.0\text{ fs}$, `rigidBonds all`), stesso cutoff ($12.0\text{ Å}$ con switching a $10.0\text{ Å}$), stesso PME, stesso termostato di Langevin a `400 K`, stesso barostato a `1.01325 bar`, stessa durata ($2.5\text{ ns}$ = `1.250.000` step) e stessa frequenza di campionamento (`dcdfreq 1000` = `2 ps`).
3. **Natura Stocastica della MD Caotica e Campionamento Multi-Seed (`Laboratory_Project_2026.pdf`, Slide 10, 15, 21; `Lecture04`, Slide 12, 15, 18–19; *Bekker et al., 2019* p. 431 & p. 436):**
   - Un sistema di $N = 21.461$ atomi accoppiati da potenziali non-lineari è intrinsecamente **caotico (effetto farfalla)**: due traiettorie che partono dalle stesse identiche coordinate atomiche $\mathbf{r}_i(0)$ ma con velocità iniziali $\mathbf{v}_i(0)$ estratte con un diverso seme del generatore di numeri casuali (**`seed` della distribuzione di Maxwell-Boltzmann**, `Lecture04` Slide 12 & 15):
     $$P(v_{i\alpha}) = \sqrt{\frac{m_i}{2\pi k_B T}} \exp\left(-\frac{m_i v_{i\alpha}^2}{2 k_B T}\right)$$
     esplorano microstati diversi dell'ensemble statistico (`Lecture04`, Slide 18) e possono mostrare tempi di innesco dell'unfolding locale leggermente diversi.
   - Per questo motivo, sia *Bekker et al. (2019)* (p. 431 & Table II p. 434) sia il Prof. Gautieri (`Laboratory_Project_2026.pdf`, Slide 15 & 21) richiedono di eseguire **più repliche indipendenti con `velocity seeds` diversi** per ciascun sistema ($n = 3$ repliche per `WT`, `Mut1`, `Mut2` e `Double`), così da poter stimare la **Media $\pm$ Deviazione Standard ($\mu \pm \sigma$)** di ogni osservabile!

---

## 2. 🎯 Cosa Cerchiamo nel Punto 4 (Gli Obiettivi Operativi)

1. **Costruzione Identica dei Domini Solvatati dei Mutanti (`04_md_simulations/`):** Generare per `Mut1`, `Mut2` e `Double` i file `.psf` e `.pdb` solvatati e ionizzati in `CHARMM36m` partendo dalle strutture 3D prodotte da `FoldX BuildModel`, verificando che non vi siano errori topologici e che $q_{\text{tot}} = 0.0000\text{ e}$.
2. **Esecuzione della Campagna Multi-Seed a `400 K` ($4\text{ sistemi} \times 3\text{ seed} = 12\text{ traiettorie}$ a `400 K`):**
   - **`WT` (`2MYO` Wild-Type):** `seed 12345` *(✅ già eseguita nel Punto 2!)* + `seed 23456` + `seed 34567`.
   - **`Mut1` (es. `K37L` in `ANK2`):** `seed 12345` + `seed 23456` + `seed 34567`.
   - **`Mut2` (es. `T103L` in `ANK4` o `H68G` in `ANK3`):** `seed 12345` + `seed 23456` + `seed 34567`.
   - **`Double` (`Mut1 + Mut2`):** `seed 12345` + `seed 23456` + `seed 34567`.
3. **Riproducibilità e Tracciabilità (`Lecture04`, Slide 27: *Production files and reproducibility*):** Archiviare in modo ordinato i file `.conf`, `.log`, `.dcd`, `.xst` e `.coor` di ciascuna replica (`<sistema>_400K_run<1..3>.dcd`) per alimentarli allo script di analisi statistica del Punto 5.

---

## 3. 🛠️ Come Lo Eseguiremo (Protocollo Programmatico Passo-Passo)

Come dettagliato nel Livello 2 di [`ROADMAP/ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md):

1. **Step 4.1 — Script Parametrico di Setup Topologico (`scripts/08_build_matched_mutants_charmm36m.tcl`):**
   - Per ciascun PDB mutato uscito da `FoldX BuildModel`:
     - Rinominazione esplicita delle istidine (`resid 40 -> HSD`, `resid 67 -> HSE`, `resid 68 -> HSD` [se non mutata in `GLY`], `resid 80 -> HSD`) e `pdbalias atom ILE CD1 CD`;
     - Filtro `protein and noh` e ricostruzione degli idrogeni con `psfgen` + `common/charmm36/top_all36_prot.rtf`;
     - Solvatazione `solvate ... -t 12` e ionizzazione `autoionize ... -sc 0.15`;
     - Scrittura automatica dei parametri PBC (`cellBasisVector1/2/3`, `cellOrigin`) nei file di configurazione NAMD `04_md_simulations/<variant>/run_seed<1..3>.conf`.
2. **Step 4.2 — Generazione Automatica dei File `.conf` Multi-Seed e Divisione del Carico HPC (`+p8`):**
   - Tutti i file `.conf` sono cloni esatti di [`02_wt_baseline/wt_400K.conf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/wt_400K.conf), in cui variano **unicamente**:
     - `structure` e `coordinates` (che puntano al `.psf`/`.pdb` della variante corrispondente);
     - `seed` (`12345` per Replica 1, `23456` per Replica 2, `34567` per Replica 3);
     - `outputName` (`04_md_simulations/<variant>/<variant>_400K_s<1..3>`).
   - **Ottimizzazione dei tempi (Divisione del lavoro Oleg + Alessia):**
     - Delle `12` traiettorie a `400 K` (`2.5 ns` ciascuna = `~4.7 ore` su `namd3 +p8` a `12.68 ns/day`), **`1` è già pronta** (`wt_400K.dcd`).
     - Le restanti **`11` traiettorie** (`~51.7 ore` totali di CPU) possono essere lanciate in batch notturno automatico e divise a metà tra il PC Linux di Oleg (`6` run $\approx 28\text{ ore}$, 2 notti) e il PC Windows di Alessia (`5` run $\approx 23\text{ ore}$, 2 notti)!

---

## 4. 💡 Cosa Abbiamo Già Trovato al Passo 4.1 (Controllo di Qualità Topologico — Gate 2 Superato!) e Controlli in Corso

Eseguendo [`scripts/08_build_matched_mutants_charmm36m.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/08_build_matched_mutants_charmm36m.tcl), abbiamo già costruito tutti i domini solvatati e ionizzati in `04_md_simulations/` e superato al **100% il Gate 2** (vedi report completo in [`04_md_simulations/matched_systems_qc_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/matched_systems_qc_report.txt) e il commento teorico in [`04_md_simulations/RISULTATI_E_COMMENTO_FASE4_SETUP.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/RISULTATI_E_COMMENTO_FASE4_SETUP.md)):

| Sistema (`04_md_simulations/`) | Atomi Totali | Molecole $\text{H}_2\text{O}$ (`TIP3P`) | Ioni $\text{Na}^+$ (`SOD`) | Ioni $\text{Cl}^-$ (`CLA`) | Carica Totale Cella | Dimensioni Box PBC ($L_x \times L_y \times L_z\text{ Å}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`WT` (`2MYO`)** | `21.461` | `6.535` | `23` | `18` | **`+0.0000 e`** | `71.840 × 58.972 × 54.840` |
| **`K37L`** (`Mut1`) | `23.068` | `7.070` | **`26`** | `20` | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`T103L`** (`Mut2`) | `23.075` | `7.070` | **`25`** | `20` | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`K37L_T103L`** (`Double`) | `23.067` | `7.068` | **`26`** | `20` | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`H68G`** *(Opz. BTC)* | `23.081` | `7.077` | **`25`** | `20` | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`K37L_H68G_T103L`** *(Opz. Tripla)* | `23.069` | `7.072` | **`26`** | `20` | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |

1. **Neutralità Elettrica e Topologia (`psfgen` & `autoionize` log):** Carica totale della cella ionizzata pari a **`0.0000 e`** per tutte le varianti (`116` termini CMAP attivi, zero coordinate mancanti; `26 Na+` compensano esattamente la carica `-6.00 e` delle varianti con `K37L`).
2. **Stabilità Termodinamica del Solutore (`grep "ENERGY:" *.log`):** Durante le run NAMD lanciate con [`04_md_simulations/run_batch_phase4.sh`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/run_batch_phase4.sh), verificheremo l'assenza di errori `Constraint failure in RATTLE algorithm`, temperatura media stabile a **`400.0 ± 3.5 K`** e volume della cella PBC equilibrato ($P = 1\text{ atm}$).
3. **Completezza dei Frame (`.dcd`):** Esattamente **`1.252` frame** salvati per ciascuna replica da `2.5 ns` (`dcdfreq 1000`), così che la finestra finale equilibrata del 30% (`frame 876..1251`, `376` frame = `0.75 ns` finali) sia rigorosamente identica tra `WT`, `Mut1`, `Mut2` e `Double`.

---

## 5. 🏆 Risultati Reali dello Screening MD `Seed 1` a `400 K` (`WT` vs `K37L`, `T103L`, `K37L_T103L`, `K37L_H68G_T103L`) e Prossimi Passi

Il 22 Novembre 2026 è stata completata con successo l'intera batteria di screening MD `Seed 1` (`12345`, `2.5 ns` ciascuna a `400 K`) tramite [`04_md_simulations/run_4_screening_400K.sh`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/run_4_screening_400K.sh). I plot e i report separati di ogni singola simulazione sono archiviati in [`05_analysis/single_runs/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs) e commentati in dettaglio in **[`05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md)**:

- **`00_WT` (`2MYO` Wild-Type):** `Q: hydrophilic-all = 0.6750 ± 0.0067`, `Q: all-all = 0.6948 ± 0.0054`, `Q_ANK2 = 0.6015`, `RMSD_BB = 6.68 ± 0.52 Å`, `Mask3AAA_RMSD = 7.27 ± 0.99 Å`, `Rgyr = 15.74 Å`, `SASA = 8838.8 Å²`, **`41` residui dello scaffold instabili** ($\langle Q_i \rangle_{400\text{K}} < 0.60$).
- **`01_K37L` (Singolo 1, `ANK2`):** `Q: hydrophilic-all = 0.7134 ± 0.0070` (**`+3.84%`**), `Q: all-all = 0.7243 ± 0.0052` (**`+2.95%`**), `Q_ANK2 = 0.6351`, `RMSD_BB = 4.17 ± 0.33 Å`, `SASA = 7941.6 Å²` (`-897 Å²`), **`25` residui instabili (`-16` residui salvati!)**.
- **`02_T103L` (Singolo 2, `ANK4`):** `Q: hydrophilic-all = 0.7047 ± 0.0066` (**`+2.97%`**), `Q: all-all = 0.7124 ± 0.0046`, `RMSD_BB = 3.77 ± 0.19 Å`, `Mask3AAA_RMSD = 3.37 ± 0.36 Å`, `SASA = 7969.2 Å²` (`-870 Å²`), **`30` residui instabili (`-11`)**.
- **`03_K37L_T103L` (Doppio `ANK2 + ANK4`):** **`Q: hydrophilic-all = 0.7178 ± 0.0075` (`+4.28%`, #1 assoluto!)**, `Q: all-all = 0.7194 ± 0.0061`, **`Rgyr = 14.76 ± 0.34 Å` (#1 assoluto per compattezza globulare)**, interfaccia `ANK1–ANK2 = 0.266` (`×2.35` vs WT), **`31` residui instabili (`-10`)**.
- **`04_K37L_H68G_T103L` (Triplo `ANK2 + ANK3 + ANK4`):** **`Q: all-all = 0.7289 ± 0.0049` (#1 assoluto, `+3.41%`)**, **`Q: hydrophilic-all = 0.7172 ± 0.0056`**, **`Q: hydrophobic-small = 0.6631 ± 0.0140` (#1 assoluto)**, **`Q_ANK2 = 0.6444` (#1 assoluto, `+4.29%`)**, **`RMSD_BB = 3.30 ± 0.42 Å` (#1 assoluto, `-51%`)**, **`Mask3AAA_RMSD = 2.83 ± 0.26 Å` (#1 assoluto, `-61%`)**, **`SASA = 7646.9 Å²` (#1 assoluto, `-1192 Å²`)**, **`23` residui instabili (`-18` residui salvati = `-44%`, #1 assoluto!)**, con arresto totale della fluttuazione del loop `ASP 60` ($\text{RMSF}_{60} = 1.50\text{ Å}$) e di `LEU 37` ($\text{RMSF}_{37} = 0.91\text{ Å}$).

### ⏭️ Prossimo Passo: Campagna di Validazione Statistica (`Seed 2` e `Seed 3` + Verifica con Pre-Equilibrazione Restrained a `300 K`)
Completato lo screening sul `Seed 1`, la campagna di validazione finale prevede:
1. L'unione dei 5 dataset `Seed 1` nei plot comparativi multi-sistema (`fig9`–`fig11`) e nella Table II di Bekker.
2. L'esecuzione delle repliche statistiche indipendenti (`Seed 2 = 23456` e `Seed 3 = 34567`) e del protocollo di verifica con **Pre-Equilibrazione Restrained a 300 K** (`constraints on` sul backbone prima dello stress a 400 K, `Lecture04` Slide 6–9).
