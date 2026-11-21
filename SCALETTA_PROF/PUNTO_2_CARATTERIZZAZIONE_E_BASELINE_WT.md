# 🔥 PUNTO 2 DELLA SCALETTA DEL PROFESSORE — "2. Characterization and baseline: Prepare the WT and run a pilot simulation. Identify candidate regions"

**Riferimento Ufficiale nella Scaletta del Corso:**  
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 9–12** (*Reference study: single-domain antibodies — Bekker, Ma & Kamiya 2019; 300 K vs 400 K vs 500 K; Lessons for the project: Local contacts can guide candidate selection*), **Slide 13** (*Aim of the project: Evaluate protein variants stability by measuring the Q-value*), **Slide 15** (*Project workflow — Step 2: Characterization and baseline. Prepare the WT and run a pilot simulation. Identify candidate regions*), **Slide 20–21** (*MD simulations and analysis: Computational budget and pilot benchmark — Benchmark first, record system size, hardware and measured simulation speed in ns/day*), **Slide 22** (*RMSF: local flexibility — Reconstruct periodic molecules and align the trajectory -> `analysisMD.tcl`*), **Slide 23** (*Q-value - Native contacts via `VMDextensions.tcl`*).

---

## 1. 📖 Introduzione al Problema e Riferimenti alle Slide di Teoria del Corso

Una volta verificata la struttura statica di `2MYO` nel Punto 1, il **Punto 2** affronta una domanda biofisica fondamentale: **come facciamo a individuare i "punti deboli" strutturali di una proteina che si innescano ad alta temperatura usando la Dinamica Molecolare (MD) classica?**

1. **Dal Modello Atomistico al Force Field Classico (`Lecture03 - Molecular Dynamics Setup.pdf`, Slide 3–30):**
   - Poiché risolvere l'equazione di Schrödinger quantistica per migliaia di atomi è computazionalmente impossibile (`Lecture03`, Slide 3–5), applichiamo l'approssimazione di Born-Oppenheimer e trattiamo i nuclei come masse puntiformi classiche collegate da molle (approccio di Meccanica Molecolare *beads-and-springs*, `Lecture01` Slide 52; `Lecture03` Slide 6–10).
   - L'energia potenziale totale del sistema $V(\mathbf{r}_1, \dots, \mathbf{r}_N)$ è definita dal **Force Field empirico** (`Lecture03`, Slide 12–28):
     $$V = \sum_{\text{bonds}} k_b (b - b_0)^2 + \sum_{\text{angles}} k_\theta (\theta - \theta_0)^2 + \sum_{\text{dihedrals}} k_\phi [1 + \cos(n\phi - \delta)] + \sum_{\text{impropers}} k_\omega (\omega - \omega_0)^2 + \sum_{i < j} \left\{ 4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right] + \frac{q_i q_j}{4\pi\varepsilon_0 r_{ij}} \right\}$$
   - Tra i force field biomolecolari (`Lecture03`, Slide 29), abbiamo scelto **CHARMM36m** (`top_all36_prot.rtf` + `par_all36m_prot.prm`), che include i termini di correzione energetica bidimensionale del backbone (**CMAP cross-terms**) calibrati per evitare il noto bias di sovra-stabilizzazione delle $\alpha$-eliche delle versioni precedenti (`CHARMM22/27`) — un aspetto decisivo per una proteina interamente a solenoide $\alpha$-elica come la Myotrophin!
2. **Modello d'Acqua Esplicito, Periodic Boundary Conditions (PBC), Cutoff, PME e Contro-ioni (`Lecture03`, Slide 32–50):**
   - In natura le proteine operano immerse in acqua salina (`Lecture03`, Slide 32). A differenza del solvente implicito continuum (`Lecture03`, Slide 36–37, che non cattura i ponti idrogeno specifici acqua–proteina che innescano il thermal unfolding), utilizziamo il modello d'acqua esplicito a 3 siti **TIP3P** (`Lecture03`, Slide 33, 35, 38).
   - Per evitare effetti di bordo superficiale con il vuoto, immergiamo il sistema in una scatola periodica replicata all'infinito nelle 3 direzioni (**Periodic Boundary Conditions, PBC**, `Lecture03` Slide 40–41).
   - Le interazioni di Lennard-Jones a corto raggio vengono troncate dolcemente con una *switching function* tra `10.0 Å` e un **cutoff di `12.0 Å`** (`Lecture03`, Slide 42–48): per la *Minimum Image Convention* (`Lecture03`, Slide 43), il padding d'acqua su ogni lato della proteina deve essere $\ge 12.0\text{ Å}$ affinché la proteina non interagisca mai con la propria immagine periodica (distanza minima tra immagini $> 24\text{ Å}$).
   - Le interazioni elettrostatiche di Coulomb a lungo raggio ($\propto 1/r$) non possono essere semplicemente troncate senza causare gravi artefatti: vengono risolte esattamente nello spazio reciproco tramite **Particle Mesh Ewald (PME)** (`Lecture03`, Slide 49), che richiede che la scatola sia **elettricamente neutra ($q_{\text{tot}} = 0$)** mediante l'aggiunta di contro-ioni ($\text{Na}^+, \text{Cl}^-$) a concentrazione fisiologica (`Lecture03`, Slide 50).
3. **Minimizzazione, Integrazione del Moto, Termostato/Barostato e Protocollo ad Alta Temperatura (`Lecture04 - Running Molecular Dynamics.pdf`, Slide 4–27; `Laboratory_Project_2026.pdf`, Slide 9–12; *Bekker et al., 2019*):**
   - Prima di integrare le equazioni di Newton $\mathbf{F}_i = m_i \mathbf{a}_i = -\nabla_i V$ (`Lecture04`, Slide 11–15), è obbligatoria una **Minimizzazione Energetica** (`Lecture04`, Slide 4–9) per rilassare gli scontri sterici di Van der Waals ($r^{-12}$) tra il solvente appena aggiunto e la proteina.
   - Per spingere il passo di integrazione temporale a **$\Delta t = 2.0\text{ fs}$** (`Lecture04`, Slide 16–17), si vincolano le lunghezze dei legami covalenti che coinvolgono gli atomi di idrogeno (`rigidBonds all`, algoritmo SHAKE/RATTLE), controllando temperatura e pressione nell'ensemble isotermo-isobarico **NPT** tramite il termostato di Langevin e il barostato a pistone (`Lecture04`, Slide 18–21).
   - **Perché simulare a `400 K` con controllo a `300 K`? (`Laboratory_Project_2026.pdf`, Slide 9–12; *Bekker et al., 2019* Fig. 2A):**
     - Misurare direttamente l'unfolding all'equilibrio vicino alla $T_m$ sperimentale ($\sim 53\text{--}56^\circ\text{C} \approx 326\text{--}329\text{ K}$ per `2MYO`) richiederebbe scale temporali dei millisecondi (`Lecture05`, Slide 2–3, 17).
     - Come dimostrato da *Bekker et al. (2019)*, a **`300 K`** tutte le proteine restano quasi congelate intorno alla struttura nativa (correlazione debole con $T_m$, $r = 0.51$), mentre a **`500 K`** l'energia termica eccessiva fa esplodere indiscriminatamente il fold ($r = 0.08$). A **`400 K` ($127^\circ\text{C}$)** si ottiene la finestra ideale ($r = 0.79\text{--}0.84$): la struttura globale regge, ma le regioni intrinsecamente più deboli perdono i loro contatti atomici nativi ($Q_i < 0.60$)!

---

## 2. 🎯 Cosa Cercavamo (Gli Obiettivi del Punto 2)

1. **Costruzione del Dominio Fisico WT:** Generare un sistema solvatato e ionizzato `CHARMM36m` fisicamente perfetto (`2myo_wb_ionized.psf` / `.pdb`) e ricavare i vettori esatti della cella periodica (`cellBasisVector1/2/3`, `cellOrigin`).
2. **Benchmark HPC e Computational Budget (`Laboratory_Project_2026.pdf`, Slide 21):** Misurare la velocità reale di simulazione in **`ns/day`** sulla nostra macchina Linux, trovare la configurazione ottimale di thread/core e pianificare il budget computazionale dell'intero progetto.
3. **Validazione Dinamica a `300 K` vs Stress Termico a `400 K`:** Dimostrare che a `300 K` il nostro modello WT della Myotrophin rimane perfettamente stabile e compatto, mentre a `400 K` subisce un unfolding parziale selettivo.
4. **Mappatura Quantitativa delle Regioni Candidate (`Laboratory_Project_2026.pdf`, Slide 12, 15, 22–23):**
   - Quale dei 4 Ankyrin Repeats (`ANK1`, `ANK2`, `ANK3`, `ANK4`) perde più contatti nativi $Q(t)$ a `400 K`?
   - Quali singoli amminoacidi dello scaffold non-mascherato scendono sotto la soglia di instabilità di Bekker (**$\langle Q_i \rangle_{400\text{K}} < 0.60$** con $N_{\text{cont}}(0) \ge 30$)?
   - Quali legami a idrogeno specifici (`hbonds_details_300K.txt` vs `400K.txt`) si rompono passando da `300 K` a `400 K`, spiegando a livello atomistico *perché* quelle regioni cedono?

---

## 3. 🛠️ Come L'Abbiamo Cercato (Metodologia, Script e Parametri Eseguiti)

1. **Generazione Topologia, Solvatazione e Ionizzazione ([`scripts/02_build_wt_charmm36m.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/02_build_wt_charmm36m.tcl)):**
   - Applicati i tautomeri verificati al Punto 1 (`resid 40 -> HSD`, `resid 67 -> HSE`, `resid 68 -> HSD`, `resid 80 -> HSD`), filtrato `protein and noh`, applicato `pdbalias atom ILE CD1 CD` e generata la topologia con `psfgen` (`common/charmm36/top_all36_prot.rtf`).
   - Solvatato con `solvate ... -t 12` (padding $12.0\text{ Å}$) e neutralizzato/ionizzato con `autoionize ... -sc 0.15` ($0.15\text{ M NaCl}$).
   - Pulito il file dei parametri dell'acqua `common/charmm36/toppar_water_ions_namd.prm` (prime 124 righe fino a `END`, rimuovendo l'istruzione CHARMM `return` non supportata da NAMD 3).
2. **Benchmark HPC Multi-Core ([`02_wt_baseline/benchmark.conf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/benchmark.conf) e [`02_wt_baseline/benchmark.log`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/benchmark.log)):**
   - Testato `NAMD 3.0.3 Linux-x86_64-multicore` sulla CPU AMD 8 Core Fisici / 16 Thread Logici (con unità vettoriali `AVX512F`/`AVX512CD`) confrontando l'esecuzione su tutti i 16 thread SMT (`namd3 +p16`) contro gli 8 core fisici puri (`namd3 +p8`).
3. **Campagna MD Notturna del Wild-Type ([`02_wt_baseline/run_overnight_wt.sh`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/run_overnight_wt.sh)):**
   - **Controllo Fisiologico a `300 K` ([`02_wt_baseline/wt_300K.conf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/wt_300K.conf)):** `1.000` step di minimizzazione + `750.000` step di produzione NPT a `300 K` e `1.01325 bar` ($\Delta t = 2\text{ fs} \rightarrow \mathbf{1.5\text{ ns}}$, salvataggio ogni `1.000` step = `2 ps` $\rightarrow$ `752` frame in `wt_300K.dcd`).
   - **Stress Termico a `400 K` ([`02_wt_baseline/wt_400K.conf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/wt_400K.conf)):** `1.000` step di minimizzazione + `1.250.000` step di produzione NPT a `400 K` e `1.01325 bar` ($\Delta t = 2\text{ fs} \rightarrow \mathbf{2.5\text{ ns}}$, `1.252` frame in `wt_400K.dcd`).
4. **Pipeline Unificata di Post-Processing ([`scripts/03_analyze_wt_baseline.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/03_analyze_wt_baseline.tcl) + [`scripts/plot_wt_baseline.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/plot_wt_baseline.py)):**
   - Unisce l'intero script ufficiale del Prof. Gautieri ([`laboratorio/analysisMD.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/analysisMD.tcl)) con il protocollo dei contatti nativi di *Bekker et al. (2019)* (`VMDextensions.tcl`):
     1. Ricomposizione PBC (`pbc wrap -compound residue -center com -centersel "protein" -all`) e allineamento sul backbone del frame 0 (`measure fit`, Slide 22).
     2. **Backbone RMSD($t$)** e **RMSF($C_\alpha$)** per residuo (`1..118`).
     3. **Hydrogen Bonds** nel tempo + **tabella completa di occupancy (%) di ogni coppia donatore-accettore** (`hbonds -type unique`, [`hbonds_details_300K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_300K.txt) e [`hbonds_details_400K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_400K.txt)).
     4. **SASA** ($r_{\text{probe}} = 1.4\text{ Å}$), **Raggio di Girazione ($R_g$)** e **% di $\alpha$-elica** (tramite STRIDE / `mol ssrecalc`).
     5. **$Q$-value nel tempo** (cutoff $4.5\text{ Å}$ sugli atomi pesanti `noh`): per le 3 classi chimiche di *Bekker et al.* (**`all-all`**, **`hydrophilic-all`**, **`hydrophobic-small`**) e risolto per ciascuno dei 4 Ankyrin Repeats (**`Q_ANK1`**, **`Q_ANK2`**, **`Q_ANK3`**, **`Q_ANK4`**).
     6. **$\langle Q_i \rangle_{400\text{K}}$ per singolo residuo (`1..118`)** mediato sugli ultimi 30% della traiettoria a `400 K` (frame `876..1251`), con esportazione di [`02_wt_baseline/unstable_residues_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/unstable_residues_report.txt) e della struttura 3D colorata nel campo Beta (`wt_400K_q_beta.pdb`).

---

## 4. 💡 Cosa Abbiamo Trovato e Considerazioni Critiche (Risultati Quantitativi)

### 4.1 Risultati del Sistema Solvatato e del Benchmark HPC (Slide 21)
- **Sistema finale (`01_structures/2myo_wb_ionized`):** **`21.461` atomi totali** (`1.815` atomi della Myotrophin con carica $-5.00\text{ e}$, `6.535` molecole d'acqua `TIP3P`, `23` ioni $\text{Na}^+$, `18` ioni $\text{Cl}^-$, carica elettrica totale esattamente **`0.0000 e`**). Dimensioni della cella PBC: **`71.840 × 58.972 × 54.840 Å`**, centro `(7.735, -8.087, 86.003)`.
- **Scoperta HPC (`+p8` vs `+p16`):**
  - Su 16 thread logici SMT (`namd3 +p16`): `0.0188 s/step` $\rightarrow$ **`9.01 ns/day`**.
  - Su 8 core fisici puri (`namd3 +p8`): `0.0133 s/step` $\rightarrow$ **`12.68 ns/day`** (**Speedup del `+40%`**!).
  - *Spiegazione tecnica per la Slide del Computational Budget:* nei kernel fortemente vettorizzati in **AVX-512** di NAMD 3 (forze di Lennard-Jones + FFT 3D del PME su 36 patch), eseguire 2 thread SMT sullo stesso core fisico genera contesa sulle pipeline FMA/AVX-512 e raddoppia l'overhead di sincronizzazione del load balancer Charm++. Fissare l'esecuzione a `+p8` ci ha fatto risparmiare ~3 ore sulla sola Fase 2 e oltre 18 ore sull'intera campagna del progetto!

### 4.2 Confronto Globale `300 K` (Controllo) vs `400 K` (Stress Termico)
I risultati delle due simulazioni reali ([`fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf) e [`fig5b_wt_analysisMD_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5b_wt_analysisMD_observables.pdf)) confermano in pieno il comportamento previsto da *Bekker et al. (2019)*:

| Osservabile (Finestra Finale Equilibrata) | `WT 300 K` (`1.5 ns`, Controllo Fisiologico) | `WT 400 K` (`2.5 ns`, Stress Termico) | Interpretazione Biofisica |
| :--- | :---: | :---: | :--- |
| **Backbone RMSD finale** | **`~2.63 Å`** (piatto e stabile) | **`~7.10 Å`** (deformazione moderata) | A `300 K` il modello CHARMM36m è perfettamente stabile; a `400 K` avviene un rilassamento/apertura locale senza denaturazione completa. |
| **Raggio di Girazione ($R_g$)** | **`~14.55 Å`** (costante) | **`~15.4 – 16.2 Å`** (espansione termica) | La proteina mantiene la sua architettura globale a solenoide ma si "gonfia" nei moduli più plastici. |
| **Contenuto di $\alpha$-Elica (STRIDE)** | **`44% – 49%`** (`52–58` residui, nativo) | **`32% – 40%`** (srotolamento parziale alle estremità delle eliche) | Alcuni capi N/C delle $\alpha$-eliche perdono i legami H $i \to i+4$ ad alta temperatura. |
| **$Q$-value: `all-all`** | **`0.725`** | **`0.693`** | Calo selettivo dei contatti nativi sotto stress termico. |
| **$Q$-value: `hydrophilic-all`** | **`0.715`** | **`0.678`** | I contatti polari/idrofilici cedono per primi (in accordo con *Bekker et al., 2019* Fig. 2C!). |
| **$Q$-value per Modulo (`ANK1..4`)** | `ANK1: 0.770`<br>`ANK2: 0.655`<br>`ANK3: 0.678`<br>`ANK4: 0.659` | **`ANK2: 0.594` (Il più fragile!)**<br>`ANK3: 0.662`<br>`ANK4: 0.675`<br>`ANK1: 0.703` | **`ANK2` (residui `35–67`) è l'anello debole primario del solenoide a `400 K`**, scendendo sotto la soglia $Q = 0.60$! |

### 4.3 Identificazione dei Residui Localmente Instabili ($\langle Q_i \rangle_{400\text{K}} < 0.60$) e Meccanismo Atomistico H-Bond
Dal report definitivo [`02_wt_baseline/unstable_residues_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/unstable_residues_report.txt) (frame `876–1251` a `400 K`, applicando i filtri di *Bekker et al.*: $\langle Q_i \rangle < 0.60$, $N_{\text{cont}}(0) \ge 30$, esclusi i 24 residui della Functional Mask e i 6 terminali `1–3, 116–118`), abbiamo trovato **28 residui instabili nello scaffold interno**, concentrati in **3 hotspot strutturali ben precisi** la cui causa atomistica è svelata dal confronto tra [`hbonds_details_300K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_300K.txt) e [`hbonds_details_400K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_400K.txt):

1. **Hotspot Dinamico 1 — L'Ingresso `TPLH` (`LYS 37`) e il Loop `GADVNA` (`ASP 60 – PRO 64`) del Repeat 2 (`ANK2`):**
   - **Residui instabili a 400 K:** **`GLY 35`** (`0.5522`), **`LYS 37`** (**`0.4845`**, $N_0 = 93$, `RMSF = 2.96 Å`), **`LEU 39`** (`0.5920`), **`GLY 46`** (`0.3842`), **`LEU 48`** (`0.4645`), **`ILE 50 – LEU 55`** (`0.4702 – 0.5599`), **`LYS 57`** (`0.5811`), **`ASP 60`** (**`0.4645`**), **`ILE 61`** (`0.4693`), **`ASN 62`** (`0.4945`), **`PRO 64`** (**`0.4406`**), **`ASP 65`** (`0.5991`).
   - **Perché cede (spiegazione tramite `hbonds_details_300K.txt` vs `400K.txt`):**
     - In `ANK3`, la posizione 3 del motivo `TPLH` ha la Treonina canonica (`THR 70`), che a `300 K` forma il ponte a idrogeno inter-repeat `THR70-Side-OG1 -> ILE61-Main-O` (`24.60%` occupancy).
     - In `ANK2`, invece della Treonina c'è **`LYS 37`** (carica $+1$, accanto a `ARG 36` carica $+1$). A `300 K`, la lunga catena di `LYS 37` si piega verso il basso formando un **ponte salino parassita con `ASP 60`** (`LYS37-NZ -> ASP60-OD2` `26.20%` + `OD1` `21.41%` = **`47.61%` occupancy totale**!).
     - A `400 K`, questa tensione fa cedere l'intera giunzione `ANK1–ANK2` (`LEU32-Main-N -> ARG36-Main-O` crolla dal **`43.62%` allo `0.08%`**), srotola l'elica interna di `ANK2` (`ALA43 -> LEU39` crolla dal **`32.58%` al `5.03%`**; `ASP44 -> HSD40` dal **`15.03%` allo `0.48%`**) e destabilizza tutto il loop `ASP 60 – PRO 64`!
2. **Hotspot Dinamico 2 — Il Loop `GADKTV` (`90–97`) e il $\beta$-Turn (`HIS 68`) del Repeat 3 (`ANK3`):**
   - **Residui instabili/deboli a 400 K:** **`HIS 68`** (`0.7139 WEAK`, `RMSF = 3.63 Å`, adiacente a `HIS 67` che scende a `0.4336`), **`LEU 86 – LEU 88`** (`0.5517 – 0.5842`), **`LYS 90`** (`0.5746`), **`GLY 91`** (`0.5797`), **`ASP 93`** (**`0.5427`**, `RMSF = 2.82 Å`), **`LYS 94`** (**`0.4553`**), **`THR 95`** (**`0.5911`**), **`LYS 97`** (`0.5063`).
   - **Perché cede:** Sia `HIS 68` ($\phi = +38.3^\circ$) sia `ASP 93` ($\phi = +35.4^\circ$) hanno angoli di Ramachandran nel quadrante sinistrorso positivo (dove dovrebbe esserci una Glicina), mentre `LYS 94` introduce una carica positiva non compensata nella tasca idrofobica `V/I` del motivo `GADVNA`.
3. **Hotspot Dinamico 3 — L'Ingresso dell'Elica C-Cap (`LEU 102 – THR 103`, `THR 108`) del Repeat 4 (`ANK4`):**
   - **Residui instabili/deboli a 400 K:** **`LEU 102`** (**`0.4942`**, $N_0 = 95$), **`THR 103`** (**`0.5075`**, $N_0 = 84$, `RMSF = 3.16 Å`), **`THR 108`** (`0.6992 WEAK`), **`ASN 110`** (`0.5509`), **`ILE 113`** (`0.5376`).
   - **Perché cede (e perché è vitale secondo *Lowe & Itzhaki, PNAS 2007*):**
     - `ANK4` è il **nucleo di folding primario** da cui parte il ripiegamento della Myotrophin.
     - A `400 K`, il legame di N-capping **`THR103-Side-OG1 -> ASP100-Main-O`** crolla dal **`15.56%` (`300 K`) a `0.00%` (`400 K`)**, e il legame elicoidale **`ALA107-Main-N -> LEU102-Main-O`** crolla dal **`18.48%` (`300 K`) allo `0.16%` (`400 K`)**, aprendo il C-cap all'acqua calda!

---

## 5. ⏭️ Come Agire Dopo nel Futuro (Passaggio al Punto 3 e Presentazione)

1. **Collegamento diretto col Punto 3 (Slide 12 di `Laboratory_Project_2026.pdf`: *"Local contacts can guide candidate selection"*):**
   - Invece di scegliere mutazioni FoldX a caso su zone già rigide, i **3 Hotspot Dinamici** emersi dalla MD a `400 K` (**`LYS 37` / loop `59–64` in `ANK2`**, **`HIS 68` / loop `93–95` in `ANK3`**, e **`THR 103` / `THR 108` in `ANK4`**) diventano i bersagli prioritari del nostro imbuto decisionale nel Punto 3!
2. **Baseline di Riferimento per il Punto 4 e 5:**
   - La traiettoria `wt_400K.dcd` (`2.5 ns`) costituisce già la **Replica 1 (`seed 12345`) del Wild-Type** per il confronto multi-seed della Fase 4 (`Q_all = 0.693`, `Q_phil = 0.678`, `Q_ANK2 = 0.594`, `28` residui instabili dello scaffold).
3. **Cosa portare nelle Slide della Presentazione del 16 Dicembre (Slide 3 & 4, ~2.5 minuti):**
   - Tabella del **Benchmark HPC (`+p8 = 12.68 ns/day`, +40% su AVX-512)** e del **Computational Budget** (Slide 21).
   - I due pannelli vettoriali già pronti in `06_presentation/`: [`fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf) (`RMSD`, `Q(t)`, `RMSF`, `<Q_i>_400K`) e [`fig5b_wt_analysisMD_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5b_wt_analysisMD_observables.pdf) (`% Alpha-Helix`, `Hbonds`, `Rgyr`, `Q_ANK1..4`), affiancati dal rendering 3D VMD di [`scripts/viz_03_wt_instability.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_03_wt_instability.tcl) colorato per `<Q_i>_400K`.
