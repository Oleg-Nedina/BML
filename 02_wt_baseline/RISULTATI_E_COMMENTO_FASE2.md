# 🔥 RISULTATI_E_COMMENTO_FASE2.md — Commento Scientifico, Risultati MD Baseline (300 K vs 400 K) e Considerazioni Teoriche della Fase 2 (`02_wt_baseline/`)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Sistema Simulato:** Myotrophin Wild-Type solvatata e ionizzata (`01_structures/2myo_wb_ionized.psf` / `.pdb`, `21.461` atomi, `CHARMM36m + TIP3P + 0.15 M NaCl`)  
**Stato della Fase 2:** ✅ **COMPLETATA E VALIDATA AL 100%** (Run notturne `wt_300K.dcd` 1.5 ns e `wt_400K.dcd` 2.5 ns concluse il 15 Novembre 2026)  
**Riferimenti alle Slide del Corso e Letteratura:**
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf): **Slide 9–13** (*Reference study: Bekker, Ma & Kamiya 2019 — 300 K vs 400 K vs 500 K; Q-value and experimental $T_m$*), **Slide 15** (*Step 2: Characterization and baseline*), **Slide 20–21** (*Computational budget and pilot benchmark*), **Slide 22** (*RMSF: local flexibility & `analysisMD.tcl`*), **Slide 23** (*Q-value - Native contacts via `VMDextensions.tcl`*).
- [`Lecture04 - Running Molecular Dynamics.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture04%20-%20Running%20Molecular%20Dynamics.pdf): **Slide 4–9** (*Energy minimization*), **Slide 10–17** (*Integrazione delle equazioni del moto e scelta del timestep $\Delta t = 2\text{ fs}$*), **Slide 18–21** (*Ensemble NPT, termostato e barostato*), **Slide 23–27** (*Equilibration, production & readiness*), **Slide 29–33** (*Analisi delle traiettorie: RMSD, RMSF, Struttura Secondaria, Hydrogen Bonds & Salt Bridges*).
- *Bekker GJ, Ma B, Kamiya N (2019).* **Protein Science**, 28(2):429–438.
- *Lowe AR, Itzhaki LS (2007).* **PNAS**, 104(8):2679–2684.

---

## 1. ⚡ Risultati del Benchmark HPC e Giustificazione del Computational Budget (`Slide 21`)

Seguendo la prescrizione della **Slide 21 di [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf)** (*"Benchmark first: Record system size, hardware and measured simulation speed in ns/day"*), abbiamo confrontato le prestazioni di `NAMD 3.0.3 Linux-x86_64-multicore` sulla nostra CPU AMD (8 Core Fisici / 16 Thread Logici SMT con estensioni vettoriali **AVX-512**, 24 GB RAM) sul sistema da `21.461` atomi ([`02_wt_baseline/benchmark.log`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/benchmark.log)):

| Configurazione NAMD 3 | Thread / Core Utilizzati | Tempo per Step ($\Delta t = 2\text{ fs}$) | Giorni per Nanosecondo (`days/ns`) | Velocità Misurata (`ns/day`) | Speedup Relativo |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`namd3 +p16` (2-way SMT)** | 16 Thread Logici | `0.0188 s/step` | `0.1110 days/ns` | **`9.01 ns/day`** | `1.00x` (Baseline) |
| **`namd3 +p8` (Core Fisici Puri)** | **8 Core Fisici** | **`0.0133 s/step`** | **`0.0788 days/ns`** | **`12.68 ns/day`** | **`1.41x (+40.7%)`** |

### Commento Architetturale HPC:
- Nei kernel *compute-bound* di NAMD 3 (calcolo vettorizzato in **AVX-512** delle forze non-bonded di Lennard-Jones entro $12\text{ Å}$ + trasformate di Fourier 3D per il **Particle Mesh Ewald** su una scatola compatta di `21.461` atomi suddivisa in 36 patch spaziali), lanciare 2 thread SMT sullo stesso core fisico provoca **contesa sulle unità vettoriali FMA/AVX-512**, inquinamento della cache L1/L2 e raddoppio dell'overhead di sincronizzazione del load balancer Charm++.
- Fissare tutte le simulazioni del progetto a **`namd3 +p8`** ha garantito uno **speedup gratuito del `+40.7%`** (`12.68 ns/day`, ovvero **`~1.89 ore per nanosecondo`**).
- **Calcolo del Budget Computazionale (`Slide 10 & 21`):**
  - Mentre *Bekker et al. (2019)* hanno simulato $7\text{ anticorpi} \times 3\text{ temperature} \times 10\text{ seed} \times 100\text{ ns} = 21\text{ }\mu\text{s}$ su un cluster GPU industriale (`Slide 10`), per il nostro progetto (`Slide 10`: *"The course project requires a smaller, explicitly justified computational budget"*) abbiamo calibrato:
    - **Fase 2 (Baseline WT):** `1.5 ns` a `300 K` (`2.8 ore`) + `2.5 ns` a `400 K` (`4.7 ore`) = **`7.5 ore` totali** (eseguite nella notte tra il 14 e il 15 Novembre, [`overnight_master.log`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/overnight_master.log)).
    - **Fase 4 (Matched Comparison a `400 K`):** `4 sistemi` (`WT`, `K37L`, `T103L`, `K37L+T103L`) $\times$ `3 velocity seeds` $\times$ `2.5 ns` = **`30.0 ns` totali** (di cui `2.5 ns` del WT Seed 1 già fatti $\rightarrow$ `27.5 ns` rimanenti $\approx 52\text{ ore}$ di CPU, ovvero ~26 ore a testa dividendole tra i PC di Oleg e Alessia).

---

## 2. 📊 Risultati Comparativi della Dinamica del Wild-Type: `300 K` (Controllo) vs `400 K` (Stress Termico)

I grafici vettoriali generati da [`scripts/plot_wt_baseline.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/plot_wt_baseline.py) ([`06_presentation/fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf) e [`06_presentation/fig5b_wt_analysisMD_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5b_wt_analysisMD_observables.pdf)) mostrano il confronto tra le due temperature su tutte le metriche di [`analysisMD.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/analysisMD.tcl) e di *Bekker et al. (2019)*:

| Metrica (`analysisMD.tcl` + `VMDextensions.tcl`) | `WT 300 K` (`1.5 ns`, Frame `526–751`) | `WT 400 K` (`2.5 ns`, Frame `876–1251`) | Variazione ($\Delta_{400\text{K}-300\text{K}}$) | Significato Biofisico (`Lecture04` Slide 29–33; `Slide 11–13`) |
| :--- | :---: | :---: | :---: | :--- |
| **Backbone RMSD (Å)** | **`2.65 ± 0.06 Å`** | **`6.98 ± 0.14 Å`** | **`+4.33 Å`** | A `300 K` la traiettoria raggiunge subito un plateau stabile intorno a `2.6 Å` (validando il setup `CHARMM36m`); a `400 K` mostra una deformazione marcata ma non un collasso casuale. |
| **Raggio di Girazione $R_g$ (Å)** | **`14.56 ± 0.07 Å`** | **`15.82 ± 0.35 Å`** | **`+1.26 Å`** | A `300 K` il volume idrodinamico è inchiodato al valore nativo; a `400 K` il solenoide si espande per l'allontanamento parziale dei loop inter-repeat. |
| **$\alpha$-Helicity STRIDE (%)** | **`46.4%` (`~55` res)** | **`35.8%` (`~42` res)** | **`-10.6%` (`-13` res)** | A `300 K` tutte le 8 eliche native sono intatte; a `400 K` circa 13 residui alle estremità delle eliche subiscono *helix fraying*. |
| **$Q$-value: `all-all`** | **`0.723 ± 0.005`** | **`0.692 ± 0.003`** | **`-0.031`** | Perdita netta di contatti atomici pesanti ($\le 4.5\text{ Å}$) ad alta temperatura. |
| **$Q$-value: `hydrophilic-all`** | **`0.714 ± 0.005`** | **`0.678 ± 0.004`** | **`-0.036`** | **Conferma la legge di Bekker et al. (2019, Fig. 2C):** la classe `hydrophilic-all` è quella che scende più in basso a `400 K`! |
| **$Q$-value: `hydrophobic-small`** | **`0.629 ± 0.011`** | **`0.621 ± 0.013`** | `-0.008` | Il core idrofobico profondo resiste più a lungo della superficie polare, impedendo la denaturazione totale a `400 K`. |
| **$Q$-value Risolto per Modulo (`ANK1..4`)** | `ANK1: 0.768`<br>`ANK2: 0.654`<br>`ANK3: 0.678`<br>`ANK4: 0.659` | **`ANK2: 0.593` (Cede!)**<br>`ANK3: 0.662`<br>`ANK4: 0.675`<br>`ANK1: 0.703` | **`ANK2: -0.061`**<br>`ANK1: -0.065`<br>`ANK3: -0.016`<br>`ANK4: +0.016` | **Il Repeat 2 (`ANK2`, residui `35–67`) è l'unico modulo a scendere sotto la soglia critica $Q = 0.60$ a `400 K` (`0.593`, fine run `0.586`)!** |

---

## 3. 🔬 Dissezione Atomistica dei "Bulloni Deboli" a `400 K`: L'Incrocio tra $\langle Q_i \rangle_{400\text{K}}$, RMSF e Occupancy H-Bond

Applicando la definizione esatta di **Bekker et al. (2019, p. 433, Tables S1–S7)** nella finestra equilibrata finale (`frame 876–1251` su `1252`, ultimi 30% della run a `400 K`):
- Escludendo i `6` residui terminali (`1–3, 116–118`), i `24` residui della Functional Mask e i residui con $< 30$ contatti nativi iniziali, troviamo esattamente **`28` residui dello scaffold interno localmente instabili ($\langle Q_i \rangle_{400\text{K}} < 0.60$)** ([`02_wt_baseline/unstable_residues_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/unstable_residues_report.txt)):
  - **`ANK2` (`13` residui instabili nello scaffold!):** `GLY 35 (0.552)`, **`LYS 37 (0.485)`**, `LEU 39 (0.592)`, `GLY 46 (0.384)`, `LEU 48 (0.465)`, `ILE 50 (0.470)`, `LEU 51 (0.480)`, `PHE 53 (0.534)`, `LEU 54 (0.531)`, `LEU 55 (0.560)`, `LYS 57 (0.581)`, **`ASP 60 (0.465)`**, `ILE 61 (0.469)`, `ASN 62 (0.495)`, **`PRO 64 (0.441)`**, `ASP 65 (0.599)`.
  - **`ANK3` (`8` residui instabili nello scaffold):** `PRO 71 (0.5995)`, `LEU 86 (0.584)`, `LEU 87 (0.577)`, `LEU 88 (0.552)`, `LYS 90 (0.575)`, `GLY 91 (0.580)`, **`ASP 93 (0.543)`**, **`LYS 94 (0.455)`**, **`THR 95 (0.591)`**, `LYS 97 (0.506)` (più **`HIS 68`** `WEAK` a `0.714` con picco di **`RMSF = 3.63 Å`**).
  - **`ANK4` (`4` residui instabili nel C-cap):** **`LEU 102 (0.494)`**, **`THR 103 (0.508)`**, `ASN 110 (0.551)`, `ILE 113 (0.538)` (più **`THR 108`** `WEAK` a `0.699`).
  - **`ANK1` (`10` residui instabili nell'N-cap):** `PHE 6 (0.563)`, `LEU 10 (0.495)`, `LEU 15 (0.576)`, `TYR 21 (0.553)`, `ASP 27 (0.449)`, `VAL 28 (0.404)`, `ASN 29 (0.574)`, `ARG 30 (0.570)`, `THR 31 (0.568)`, `LEU 32 (0.397)`, `GLU 33 (0.593)`.

### La Prova Regina dalla Tabella degli Hydrogen Bonds (`hbonds_details_300K.txt` vs `hbonds_details_400K.txt`):
Confrontando l'occupancy (%) dei legami a idrogeno intramolecolari tra `300 K` ([`hbonds_details_300K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_300K.txt)) e `400 K` ([`hbonds_details_400K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_400K.txt)), abbiamo ricostruito la **catena cinematica esatta del cedimento termico**:

| Legame a Idrogeno / Ponte Salino | Modulo / Funzione Strutturale | Occupancy a `300 K` | Occupancy a `400 K` | Cosa Rivelano Questi Numeri sul Meccanismo di Unfolding |
| :--- | :--- | :---: | :---: | :--- |
| **`LYS37-Side-NZ -> ASP60-Side-OD1/OD2`** | Ponte salino parassita intra-`ANK2` | **`47.61%`** (`26.20% + 21.41%`) | `25.88%` (`17.73% + 8.15%`) | Invece di stabilizzare il motivo `TPLH` (`37–40`), `LYS 37` tira verso di sé `ASP 60` deformando il loop `GADVNA` di `ANK2`! |
| **`LEU32-Main-N -> ARG36-Main-O`** | Cucitura inter-repeat `ANK1–ANK2` | **`43.62%`** | **`0.08%` (Rotto!)** | A `400 K` la giunzione tra `ANK1` e `ANK2` (proprio accanto a `LYS 37`) si stacca completamente! |
| **`ALA43-Main-N -> LEU39-Main-O`** | H-bond $i \to i+4$ Elica interna `ANK2` (`TPLH`) | **`32.58%`** | **`5.03%` (Crollato!)** | Priva di un N-cap stabile in posizione 37, l'elica interna di `ANK2` (`HELIX 3`) si srotola a `400 K`. |
| **`ASP44-Main-N -> HSD40-Main-O`** | H-bond $i \to i+4$ Elica interna `ANK2` | **`15.03%`** | **`0.48%` (Rotto!)** | Lo srotolamento di `HELIX 3` arriva a minacciare il residuo funzionale `ASP 44`! |
| **`THR70-Side-OG1 -> ILE61-Main-O`** | Ponte canonico `TPLH` `ANK2–ANK3` | **`24.60%`** | **`0.08%` (Rotto!)** | Il cedimento del loop `ASP 60 – ILE 61 – PRO 64` di `ANK2` fa saltare l'aggancio con `THR 70` di `ANK3`. |
| **`THR103-Side-OG1 -> ASP100-Main-O`** | N-capping dell'elica C-cap `ANK4` | **`15.56%`** | **`0.00%` (Rotto!)** | A `400 K` l'ingresso del C-cap (`THR 103`, $\langle Q_{103} \rangle = 0.508$) perde il contatto di capping con `ASP 100`. |
| **`ALA107-Main-N -> LEU102-Main-O`** | H-bond elicoidale del C-cap `ANK4` | **`18.48%`** | **`0.16%` (Rotto!)** | L'inizio dell'elica C-cap (`LEU 102 – THR 103`) si apre all'acqua calda. |

---

## 4. 🧠 Considerazioni Teoriche: Perché RMSF da Solo Non Basta (`Slide 22`) e Come Questa Fase Guida la Fase 3

Nella **Slide 22 di [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf)**, il Prof. Gautieri scrive un avvertimento fondamentale in rosso:
> *"NOTE: flexibility alone does not establish local instability or identify a beneficial mutation."*

I nostri dati di Fase 2 dimostrano in modo esemplare la verità di questa affermazione:
1. Se avessimo guardato **solo i picchi di RMSF a `400 K`** ([`fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf) quadrante C), i valori più alti in assoluto si trovano sulle code terminali (`MET 1`, `GLN 118`, `RMSF > 4.5 Å`) e sulle punte dei $\beta$-hairpin funzionali (`HIS 67 – ILE 69`, `RMSF ~ 3.6 Å`), che però devono restare flessibili o intoccabili per legare l'Actin-Capping Protein!
2. Guardando invece la **perdita di contatti nativi $\langle Q_i \rangle_{400\text{K}} < 0.60$ unita alla rottura dei legami a idrogeno elicoidali** (`Slide 12 & 23`), abbiamo scoperto che i veri cedimenti strutturali interni avvengono all'**ingresso del motivo `TPLH` di `ANK2` (`LYS 37`, $\langle Q_{37} \rangle = 0.485$)** e all'**ingresso dell'elica C-cap di `ANK4` (`LEU 102 – THR 103`, $\langle Q_{103} \rangle = 0.508$)**, che hanno un RMSF moderato (`~3.0 Å`) ma perdono oltre il **50% dei loro contatti atomici nativi**!
3. Questo spiega perché le mutazioni **`K37L`** (in `ANK2`) e **`T103L`** (in `ANK4`) selezionate nella Fase 3 non sono semplici "punteggi energetici di un software", ma rispondono esattamente alla diagnosi fisica della Dinamica Molecolare del Punto 2!
