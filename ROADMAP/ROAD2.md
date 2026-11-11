# 🗺️ ROAD2.md — Fase 2: Benchmark HPC, Computational Budget e Dinamica Molecolare Baseline del Wild-Type (`2MYO`)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Target:** Myotrophin Wild-Type (`2MYO`), sistema solvatato e ionizzato `01_structures/2myo_wb_ionized.psf` / `.pdb` (`21.461` atomi)  
**Cartella di lavoro di questa fase:** `02_wt_baseline/`  
**Regole di Riferimento:** Tutte le operazioni seguono tassativamente quanto stabilito in [`AGENTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/AGENTS.md) (*Zero Ipotesi, Citazione Obbligatoria della Documentazione, Esecuzione da parte dell'Utente, Note Foto per le Slide*).

---

## 📚 1. Fonti e Documentazione di Riferimento per la Fase 2

Ogni parametro e comando utilizzato in questa fase deriva direttamente da:

1. **Linee Guida Ufficiali del Progetto:**  
   [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 11–13** (*Why elevated temperature? Finding the right temperature — Bekker et al. 300 K vs 400 K vs 500 K*), **Slide 15** (*Workflow Step 2: Benchmark & WT Simulation*), **Slide 21** (*Computational Budget & Replicas via Velocity Seeds*), **Slide 22** (*Trajectory Alignment & RMSF vs Local Instability*).
2. **Slide delle Lezioni del Corso (Prof. A. Gautieri):**
   - [`Lecture04 - Running Molecular Dynamics.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture04%20-%20Running%20Molecular%20Dynamics.pdf) — **Slide 14–25** (Integrazione di Verlet, timestep $\Delta t = 2.0\text{ fs}$ con `rigidBonds all`, cutoff $12.0\text{ Å}$ con `switching on` a $10.0\text{ Å}$, PME per l'elettrostatica, termostato di Langevin e barostato di Nosé-Hoover Langevin Piston), **Slide 26–30** (Minimizzazione, Equilibrazione, Produzione e lettura del `.log`), **Slide 31–34** (Analisi di RMSD e RMSF).
3. **File di Configurazione e Script dei Tutorial del Corso:**
   - [`Tutorials/namd-tutorial-files/1-3-box/ubq_wb_eq.conf`](file:///home/olly/UNI/BioMolEng2026/Tutorials/namd-tutorial-files/1-3-box/ubq_wb_eq.conf) (Template ufficiale NAMD del corso per proteine in Water Box esplicito con PME e Langevin Piston).
   - [`Tutorials/nativecontacts/native_contacts.tcl`](file:///home/olly/UNI/BioMolEng2026/Tutorials/nativecontacts/native_contacts.tcl) (Calcolo della frazione di contatti nativi $Q$ con `prepareNativeContacts` e `measureNativeContacts`).
   - [`Tutorials/rmsf.tcl`](file:///home/olly/UNI/BioMolEng2026/Tutorials/rmsf.tcl) (Calcolo dell'RMSF dei carboni $\alpha$).
4. **Articolo Scientifico di Riferimento:**
   - **Bekker GJ, Ma B, Kamiya N (2019).** *Thermal stability of single-domain antibodies estimated by molecular dynamics simulations.* **Protein Science**, 28(2):429–438 (Scelta di $T = 400\text{ K}$ rispetto a $300\text{ K}$, definizione di residuo localmente instabile $\langle Q_i \rangle < 0.6$ con almeno 30 contatti nativi sugli atomi pesanti entro $4.5\text{ Å}$).

---

## 🎯 2. Obiettivi Operativi della Fase 2 (`ROAD2` — Tradotti per HPC)

In questa fase passiamo dal pre-processing all'**esecuzione del solutore parallelo NAMD 3** e al **post-processing delle traiettorie**:

1. **Passo 2.1 — Benchmark HPC (`ns/day`) sul tuo processore a 16 thread:**  
   Lanciamo un test breve (`2.500` step = $5\text{ ps}$) del sistema da `21.461` atomi con `namd3 +p16` per misurare la velocità reale di calcolo in **nanosecondi simulati al giorno (`ns/day`)** e lo speedup sui core.
2. **Passo 2.2 — Definizione del Computational Budget (Slide 21 del Progetto):**  
   In base ai `ns/day` misurati al Passo 2.1, fissiamo la durata esatta $t_{\text{sim}}$ delle nostre traiettorie e il numero di repliche (`seeds`) in modo che l'intero progetto sia fattibile nei tempi previsti ma statisticamente solido.
3. **Passo 2.3 — Simulazione MD Baseline del Wild-Type a `300 K` (Controllo Fisiologico) e a `400 K` (Stress Termico):**  
   Eseguiamo:
   - **Minimizzazione energetica (`minimize 1000`):** rilassa le sovrapposizioni iniziali tra gli atomi d'acqua appena aggiunti e la proteina (gradiente coniugato).
   - **Dinamica Molecolare in NPT ($P = 1\text{ atm}$):** una run a **$300\text{ K}$** ($27^\circ\text{C}$, dove la proteina deve restare stabile con $Q \approx 0.90\text{--}0.95$) e una run ad alta temperatura a **$400\text{ K}$** ($127^\circ\text{C}$, dove l'energia cinetica più alta mette sotto sforzo i legami deboli e fa cedere i "bulloni" più fragili della proteina!).
4. **Passo 2.4 — Post-Processing della Traiettoria (Allineamento PBC, RMSD, RMSF, $Q$-value Globale, per Repeat `ANK1–4` e per singolo Residuo $Q_i$):**  
   Identifichiamo matematicamente quali dei 94 amminoacidi "non-read-only" perdono i loro contatti nativi a $400\text{ K}$ ($\langle Q_i \rangle < 0.6$). Quelli saranno i candidati che muteremo nella Fase 3!

---

## ⚠️ Nota Tecnica Importante sul File `toppar_water_ions.str` (Per far girare NAMD con `CHARMM36m`)

Nel Passo 1.5 abbiamo generato il file `.psf` con la topologia moderna **`CHARMM36m`** (`common/charmm36/top_all36_prot.rtf`).  
Quando NAMD avvia la simulazione, deve leggere i parametri delle forze da:
1. `common/charmm36/par_all36m_prot.prm` (forze della proteina `CHARMM36m`)
2. `common/charmm36/toppar_water_ions.str` (forze dell'acqua `TIP3P` e degli ioni `SOD`/`CLA`)

Ispezionando il file `Tutorials/charmm36/toppar_water_ions.str` fornito nei file del corso (che abbiamo copiato in `common/charmm36/toppar_water_ions.str`), si nota che **dopo la riga 124 (`END`)** sono rimaste delle righe pensate per i lipidi (`OCL`, `OBL`) e l'istruzione finale `return` (riga 194) del linguaggio CHARMM che **NAMD non riconosce**!  
Per evitare che NAMD si fermi con l'errore `UNKNOWN PARAMETER: return` alla riga 194, basta tenere solo le prime 124 righe di `common/charmm36/toppar_water_ions.str` (che contengono tutti i parametri dell'acqua `TIP3P` e degli ioni `SOD`/`CLA` fino a `END`).

### 💻 Comando Preliminare 2.0 (Pulizia delle righe lipidi/`return` dopo la riga 124 di `toppar_water_ions.str`):

Dalla cartella `BML/`, esegui questo comando per troncare `common/charmm36/toppar_water_ions.str` alla riga 124 (`END`):

```bash
head -n 124 common/charmm36/toppar_water_ions.str > common/charmm36/toppar_water_ions_namd.prm
```

Useremo `common/charmm36/toppar_water_ions_namd.prm` nei file `.conf` di NAMD: è 100% compatibile con NAMD 3 senza alcun warning!

---

## ⚡ PASSO 2.1 — Benchmark HPC di NAMD 3 (`ns/day`)

Come spiegato nella **Slide 15 e Slide 21 di [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf)**, prima di decidere quanto lunghe fare le simulazioni dobbiamo **misurare la velocità reale di NAMD (`ns/day`)** sul nostro sistema da `21.461` atomi.

### Spiegazione dei Parametri del Solutore NAMD (Da `Tutorials/namd-tutorial-files/1-3-box/ubq_wb_eq.conf` e `Lecture04`):
- **`timestep 2.0` + `rigidBonds all`:** Ogni passo d'integrazione avanza il tempo di **$2\text{ femtosecondi}$ ($2 \times 10^{-15}\text{ s}$)**. Per poter usare $2\text{ fs}$ senza instabilità numerica, blocchiamo la vibrazione ad altissima frequenza delle molle Idrogeno–Atomo pesante (`rigidBonds all`, algoritmo SHAKE/RATTLE). Con $\Delta t = 2\text{ fs}$, **$500.000\text{ step} = 1\text{ nanosecondo (ns)}$**.
- **`cutoff 12.0`, `switchdist 10.0`, `pairlistdist 14.0`:** Le forze di Van der Waals vengono azzerate dolcemente tra $10\text{ Å}$ e $12\text{ Å}$; la *neighbor list* delle coppie di particelle vicine viene calcolata entro $14\text{ Å}$.
- **`PME yes`, `PMEGridSpacing 1.0`:** L'elettrostatica di Coulomb a lungo raggio viene risolta esattamente su una griglia 3D con passo $1\text{ Å}$ tramite FFT (*Particle Mesh Ewald*), valutata ogni 2 step (`fullElectFrequency 2`, schema *r-RESPA multi-timestep*).
- **`langevin on` + `langevinPiston on`:** Termostato stocastico di Langevin ($\gamma = 1\text{ ps}^{-1}$) e barostato a pistone ($P = 1.01325\text{ bar} = 1\text{ atm}$) per simulare nell'ensemble **NPT** (numero di particelle, pressione e temperatura costanti).
- **`cellBasisVector1/2/3` e `cellOrigin`:** I valori esatti che avete appena misurato al Passo 1.6 (`71.840`, `58.972`, `54.840` e centro `7.735, -8.087, 86.003`)!

### 💻 Comando 2.1 — Lanciare il Benchmark di NAMD 3 sui tuoi 16 thread

Una volta creato il file `02_wt_baseline/benchmark.conf`, lancialo dal terminale (dalla cartella `BML/`) con `+p16` (usa tutti i 16 thread della tua CPU) e leggi la velocità in `days/ns` e `ns/day`:

```bash
# 1. Lancia il benchmark (100 step di minimizzazione + 2500 step = 5 ps di MD)
namd3 +p16 02_wt_baseline/benchmark.conf > 02_wt_baseline/benchmark.log

# 2. Estrai le righe di Benchmark time e calcola i ns/day!
grep -E "Info: Benchmark time:|TIMING:" 02_wt_baseline/benchmark.log
```

> **👉 Come leggere l'output del Benchmark (Mentalità HPC):**  
> NAMD stamperà 3 righe `Info: Benchmark time: 16 CPUs ... s/step ... days/ns`.  
> Dividendo $1 / (\text{days/ns})$ otteniamo esattamente quanti **nanosecondi al giorno (`ns/day`)** macina il tuo PC! *(Se volete fare una mini-tabella di **Strong Scaling HPC** da mettere nella Slide 3 — che al prof piacerà tantissimo visto che venite da HPC! — potete provare a lanciarlo con `+p4`, `+p8` e `+p16` per vedere come scala sui core fisici vs HyperThreading!).*

---

## 📊 PASSO 2.2 — Calcolo del Computational Budget (Slide 21)

Appena vediamo i `ns/day` del Passo 2.1, compileremo insieme questa tabella (richiesta esplicitamente nella **Slide 21 di `Laboratory_Project_2026.pdf`**):

$$\text{Tempo Totale Simulato} = N_{\text{sistemi}} \times N_{\text{temperature}} \times N_{\text{repliche (seeds)}} \times t_{\text{traiettoria}}$$

| Voce del Budget | Valore Pianificato | Giustificazione Scientifica |
| :--- | :--- | :--- |
| **Hardware & Parallelismo** | CPU 16 threads (`AVX-512`), `NAMD 3.0.3 Multicore` | Velocità misurata dal benchmark: **`X.XX ns/day`** |
| **Dimensione Sistema** | `21.461` atomi (`6.535` $\text{H}_2\text{O}$ TIP3P + `23` $\text{Na}^+$ + `18` $\text{Cl}^-$) | Solvente esplicito PME con padding $12\text{ Å}$ e $0.15\text{ M NaCl}$ |
| **Fase 2 (Baseline WT)** | $1 \times 300\text{ K}$ + $1 \times 400\text{ K}$ | Identificazione dei residui localmente instabili ($\langle Q_i \rangle < 0.6$) |
| **Fase 4 (Matched Comparison)** | **4 sistemi** (`WT`, `Mut1`, `Mut2`, `Double`) $\times$ **3 repliche (`seeds`)** a $400\text{ K}$ | Calcolo rigoroso di Media $\pm$ Deviazione Standard ($\mu \pm \sigma$) |

---

## 🔥 PASSO 2.3 — Simulazioni MD del Wild-Type: Controllo a `300 K` vs Stress Termico a `400 K`

### Perché simuliamo il WT sia a `300 K` che a `400 K`? (In parole semplici)
- A **`300 K` ($27^\circ\text{C}$)**: è la temperatura ambiente/fisiologica. Ci serve come **controllo di sanità (sanity check)**: la proteina deve rimanere ferma nella sua forma nativa ($\text{RMSD} \approx 1.5\text{--}2.5\text{ Å}$ e frazione di contatti nativi $Q > 0.85\text{--}0.90$). Se una proteina si smontasse già a $300\text{ K}$, vorrebbe dire che abbiamo sbagliato il setup!
- A **`400 K` ($127^\circ\text{C}$)**: come dimostrato da **Bekker et al. (2019)** (e mostrato nella **Slide 11–13 di `Laboratory_Project_2026.pdf`**), a $300\text{ K}$ servirebbero millisecondi per vedere chi cede, mentre a $500\text{ K}$ esplode tutto subito. A **$400\text{ K}$**, invece, l'agitazione termica accelera i movimenti quel tanto che basta per far **allentare i "bulloni" (amminoacidi) più deboli** in pochi nanosecondi, lasciando però intatto il telaio generale!

### Protocollo di Esecuzione (Per `wt_300K.conf` e `wt_400K.conf`):
1. **Minimizzazione (`minimize 1000`):** 1000 passi di *Conjugate Gradient* a $T = 0\text{ K}$ per eliminare qualsiasi collisione tra gli atomi d'acqua aggiunti da `solvate` e la proteina.
2. **Inizializzazione delle velocità (`reinitvels $temperature`):** Assegna a ogni atomo una velocità casuale estratta dalla distribuzione di Maxwell-Boltzmann alla temperatura target ($300\text{ K}$ o $400\text{ K}$) usando un `seed` fissato e documentato (`seed 10101`).
3. **Produzione NPT (`run ...`):** Integrazione delle equazioni del moto di Newton salvando i fotogrammi del "film" nel file binario `.dcd` ogni $2\text{ ps}$ (`dcdfreq 1000`) o $5\text{ ps}$ (`dcdfreq 2500`).

---

## 📈 PASSO 2.4 — Post-Processing della Traiettoria WT e Caccia ai "Bulloni Deboli" ($Q_i$ e RMSF)

Una volta terminate le due simulazioni del WT (`300 K` e `400 K`), lo script di analisi automatica (`scripts/03_analyze_wt_baseline.tcl`) eseguirà i 4 passaggi prescritti dalla **Slide 22 di `Laboratory_Project_2026.pdf`** e dal paper di **Bekker et al. (2019)**:

1. **Ricostruzione PBC e Allineamento Roto-Traslazionale (Obbligatorio!):**
   - Durante la simulazione, la proteina fluttua e ruota liberamente dentro la scatola d'acqua. Prima di misurare di quanto si sono mossi i suoi atomi, usiamo `pbc wrap` e sovrapponiamo (`measure fit`) tutti i fotogrammi del film sul fotogramma 0, così eliminiamo la rotazione/traslazione globale della proteina e misuriamo **solo la deformazione interna**!
2. **Calcolo di RMSD($t$) e RMSF per residuo (`02_wt_baseline/rmsf_wt_300K.dat` vs `rmsf_wt_400K.dat`):**
   - **RMSD (Root Mean Square Deviation):** è una curva nel tempo che dice *"di quanti Å, in media, tutta la proteina si è deformata rispetto all'istante zero"*.
   - **RMSF (Root Mean Square Fluctuation):** è un grafico con i **118 amminoacidi sull'asse X** e **quanto oscilla ciascun amminoacido sull'asse Y**. Dove vedremo i picchi alti a $400\text{ K}$ rispetto a $300\text{ K}$, lì la struttura sta vibrando forte!
3. **Calcolo del $Q$-value Globale e per singolo Modulo Ankyrin (`ANK1`, `ANK2`, `ANK3`, `ANK4`):**
   - Il **$Q$-value (Fraction of Native Contacts)** è il numero chiave di tutto il progetto:
     - All'istante 0 contiamo tutte le coppie di atomi che si toccano (distanza $\le 4.5\text{ Å}$ sugli atomi pesanti, oppure $\le 7.0\text{ Å}$ sui $C_\alpha$). All'inizio $Q(0) = 1.00$ ($100\%$ dei contatti intatti).
     - Man mano che la proteina si scalda a $400\text{ K}$, alcuni contatti si staccano e $Q(t)$ scende (es. a $0.75$, cioè il $75\%$ dei contatti originali è ancora intatto).
     - Calcolandolo separatamente per **`ANK1`**, **`ANK2`**, **`ANK3`** e **`ANK4`**, vedremo **quale dei 4 moduli della molla è il più fragile al calore**!
4. **Calcolo del $Q$-value per singolo amminoacido ($\langle Q_i \rangle$, Tabella S1 di Bekker et al.):**
   - Per ognuno dei 118 amminoacidi, calcoliamo la percentuale media dei suoi contatti nativi rimasti intatti nella parte finale della simulazione a $400\text{ K}$.
   - Secondo la definizione di **Bekker et al. (2019, pag. 432)**, ogni amminoacido (con almeno 30 contatti nativi iniziali e fuori dalla *Functional Mask*) che scende sotto **$\langle Q_i \rangle < 0.60$** (cioè perde più del $40\%$ dei suoi agganci!) è classificato ufficialmente come **"Unstable Residue" (Bullone Instabile)**!

---

## 📸 MOMENTO FOTO PER LE SLIDE (Fase 2 — Da salvare in `06_presentation/`)

Dalla Fase 2 ricaveremo **2 figure chiave per la Slide 4 (Caratterizzazione Baseline del Wild-Type)**:

1. **Figura 4 (`06_presentation/fig4_wt_300K_vs_400K_superposition.png`) — Sovrapposizione 3D in VMD del WT a 300 K vs 400 K:**
   - Uno script grafico (`scripts/viz_03_wt_instability.tcl`) caricherà l'ultimo fotogramma a $300\text{ K}$ (in Blu trasparente) sovrapposto all'ultimo fotogramma a $400\text{ K}$ colorato in base alla perdita di contatti $\langle Q_i \rangle$ (dal Blu = stabile al Rosso = instabile, usando il campo `Beta` di VMD come in `Tutorials/beta.tcl`!).
   - In un solo colpo d'occhio si vedrà quali eliche/loop si sono mossi a $400\text{ K}$!
2. **Figura 5 (`06_presentation/fig5_wt_baseline_plots.png`) — Grafici Python Multi-Pannello (`RMSD`, `RMSF 300K vs 400K`, `Q(t)` per Repeat e `Q_i` per residuo):**
   - Uno script Python (`scripts/plot_wt_baseline.py`) genererà automaticamente il pannello ad alta risoluzione (300 DPI) con le fasce dei 4 Ankyrin Repeats e la linea di soglia $Q_i = 0.6$!

---

## ✅ Checklist Operativa della Fase 2 (`ROAD2`)

- [ ] **Passo 2.0:** Creare `common/charmm36/toppar_water_ions_namd.prm` (privo delle righe lipidi/`return` dopo la riga 124).
- [ ] **Passo 2.1:** Eseguire il benchmark NAMD (`02_wt_baseline/benchmark.conf`) con `+p16` (ed eventuale test di scaling `+p4`, `+p8`, `+p16`) e annotare i **`ns/day`**.
- [ ] **Passo 2.2:** Fissare insieme il **Computational Budget** in base ai `ns/day` misurati.
- [ ] **Passo 2.3:** Lanciare le simulazioni baseline del WT a **`300 K`** e **`400 K`**.
- [ ] **Passo 2.4:** Eseguire lo script di analisi (`RMSD`, `RMSF`, $Q$ globale, $Q$ per Repeat, $Q_i$ per residuo), generare le Figure 4 e 5 e leggere la classifica dei residui instabili ($\langle Q_i \rangle < 0.6$) da passare a FoldX/Rosetta/FireProt nella Fase 3!
