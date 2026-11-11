# 📜 AGENTS.md — Regole di Progetto, Rigore Scientifico e Buone Norme (BioMolEng 2026)

Questo documento definisce le **regole tassative e inderogabili** per lo sviluppo del progetto di **Biomolecular Modelling Laboratory (2026)** del **Team #10 (Oleg Nedina & Alessia Rigoni)** sulla proteina **Myotrophin (PDB ID: `2MYO`)**.  
Qualsiasi assistente AI (Antigravity / LLM) e ciascun membro del team deve attenersi rigorosamente a queste direttive in ogni fase del lavoro.

---

## 🚨 1. Regole Fondamentali di Interazione e Operatività

1. **Ruolo Didattico e Scientifico (L'Utente Esegue i Comandi):**
   - L'assistente AI agisce come **mentore scientifico e assistente di programmazione**, **MAI come esecutore autonomo**.
   - **I comandi da terminale, gli script VMD/TCL, le esecuzioni di NAMD3, FoldX e Python devono essere eseguiti dall'utente** (Oleg su Linux e Alessia su Windows). L'assistente fornisce i blocchi di comandi spiegati riga per riga, attendendo l'output dell'utente per commentarlo e validarlo insieme.
2. **Organizzazione Modulare nella Cartella `ROADMAP/` (`ROAD1.md`, `ROAD2.md`, ...):**
   - Il progetto avanza tramite documenti operativi incrementali salvati nella cartella [`ROADMAP/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP):
     - `ROAD1.md`: Fase 1 — Analisi Biofisica e Strutturale di `2MYO`, Protonazione Istidine, Functional Masking (`3AAA`), Consenso Ankyrin e Setup Topologico/Solvatazione (`psfgen`, `solvate`, `autoionize`).
     - `ROAD2.md`: Fase 2 — Benchmark Computazionale (`ns/day`), Computational Budget, Simulazioni MD Baseline del Wild-Type (300 K vs 400 K) e Mappatura delle Vulnerabilità Dinamiche ($\langle Q_i \rangle$, RMSF, $Q$ per Ankyrin Repeat).
     - `ROAD3.md`: Fase 3 — Multi-Tool Consensus Screening (*FoldX + FireProt 2.0 + Rosetta + ThermoMPNN/ProteinMPNN + DynaMut2 + Ankyrin Consensus*) e Costruzione 3D dei Mutanti (`BuildModel`).
     - `ROAD4.md`: Fase 4 — Simulazioni MD Comparative (*Matched Comparison*) Multi-Seed per WT, Mutante 1, Mutante 2 e Doppio Mutante.
     - `ROAD5.md`: Fase 5 — Analisi Statistica Avanzata ($Q$-value globale, chimico, per-repeat, ciclo epistatico del doppio mutante, $\Delta C_{ij}$, H-bond/Salt Bridges, elicità) e Preparazione Presentazione Finale (16 Dicembre).
   - Ogni file `ROAD<N>.md` viene generato/aggiornato prima di iniziare la fase corrispondente e deve contenere: **Obiettivi**, **Riferimenti alla Documentazione/Lezioni**, **Comandi Passo-Passo (Linux + Windows)**, **Cosa Controllare nell'Output (Quality Control)** e **Criteri di Superamento della Fase**.
3. **Zero Ipotesi ("No Guessing / No Hallucinations"):**
   - È **severamente vietato** tirare a indovinare numeri di residuo, stati di protonazione, parametri di force field, dimensioni del box o risultati energetici.
   - Ogni affermazione sulla struttura `2MYO` o `3AAA` deve derivare dalla lettura diretta dei record ufficiali del file PDB (`HEADER`, `REMARK`, `DBREF`, `SEQRES`, `HELIX`, `ATOM`) o da misure quantitative effettuate in VMD.
4. **Citazione Obbligatoria delle Fonti e della Documentazione:**
   - Ogni procedura, comando o parametro deve citare esplicitamente la propria fonte:
     - **Lezioni del Corso:** `Lecture01`–`Lecture05` del Prof. Alfonso Gautieri e [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) (indicando il numero esatto della slide).
     - **Tutorial Ufficiali del Corso:** `/home/olly/UNI/BioMolEng2026/Tutorials/` (`charmm36/`, `nativecontacts/`, `foldx_test/`, `namd-tutorial-files/`, `dhfr/`, `rmsf.tcl`).
     - **Documentazione Tecnica Ufficiale:** Manuale Utente NAMD 3.0, VMD User's Guide (`psfgen`, `solvate`, `autoionize`, `pbctools`), Documentazione `VMDextensions.tcl` (T. Giorgino), Manuale FoldX 5 Suite.
     - **Letteratura Scientifica Peer-Reviewed:**
       - *Bekker, Ma & Kamiya (2019)*, Protein Science 28:429–438 (protocollo MD ad alta temperatura e $Q$-value).
       - *Yang, Nanduri, Sen & Qin (1998)*, Structure 6:619–626 (struttura NMR di `2MYO`).
       - *Takeda et al. (2010)*, PLoS Biology 8:e1000416 (complesso cristallografico V-1 / Actin-Capping Protein, PDB `3AAA`).
       - *Mosavi, Minor & Peng (2002)*, PNAS 99:16029–16034 & *Lowe & Itzhaki (2007)*, PNAS 104:2679–2684 (consenso strutturale e meccanismo di folding degli Ankyrin Repeats di Myotrophin).

---

## 🧪 2. Buone Norme del Settore (Computational Biophysics & Protein Engineering Best Practices)

### A. Preparazione Strutturale e Chimico-Fisica (*"Garbage In, Garbage Out"*)
1. **Distinzione Metodo Sperimentale (NMR vs X-Ray):**
   - `2MYO` è una struttura **Solution NMR** (*minimized average structure*, 1 modello, determinata a $T = 297\text{ K}$, $\text{pH} = 6.2$, `REMARK 210`).
   - Contiene già atomi di idrogeno con nomenclatura NMR (`H`, `H1`, `H2`, `H3`, `HB2`, `HB3`...) e non contiene acqua cristallografica. Prima di usare `psfgen`, è obbligatorio filtrare **`protein and noh`** (solo atomi pesanti della proteina) affinché `psfgen` ricostruisca gli idrogeni secondo la topologia standard CHARMM (`Lecture03`, Slide 20–24; `Tutorials/nativecontacts/psfgen.tcl`).
2. **Verifica della Completezza e della Numerazione (`DBREF` & `SEQRES`):**
   - Verificare sempre l'assenza di gap interni nella catena polipeptidica (`getFasta` in `VMDextensions.tcl` non deve mostrare interruzioni `/`).
   - Tenere sempre traccia dello shift di numerazione tra PDB (`1..118`, dove `MET 1` è presente) e UniProt (`P62775`, dove `MET 1` iniziale è scisso, `DBREF` riga 163 di `2MYO.pdb`).
3. **Assegnazione Razionale dello Stato di Protonazione (`HIS`, `CYS`, `ASP`, `GLU`, `LYS`, `ARG`):**
   - Nel force field CHARMM, il nome residuo `HIS` non è ammesso: deve essere convertito in **`HSD`** (neutra, protone su $N_{\delta 1}$), **`HSE`** (neutra, protone su $N_{\epsilon 2}$) o **`HSP`** (protonata carica $+1$ su entrambi gli azoti) (`Lecture03`, Slide 22).
   - Non applicare mai `pdbalias residue HIS HSE` alla cieca senza aver prima ispezionato in VMD l'intorno chimico (donatori/accettori di H-bond entro $3.5\text{ Å}$) di tutte le istidine della proteina (**`HIS 40`**, **`HIS 67`**, **`HIS 68`**, **`HIS 80`** in `2MYO`).
   - Verificare le distanze $S_\gamma \text{--} S_\gamma$ tra tutte le cisteine (**`CYS 2`**, **`CYS 45`**, **`CYS 83`** in `2MYO`) per confermare l'assenza di ponti disolfuro (`DISU`).
4. **Vincoli Funzionali (*Functional Masking*):**
   - In ingegneria proteica, aumentare la termostabilità ($T_m$) non deve mai distruggere la funzione biologica (`Laboratory_Project_2026.pdf`, Slide 9 & Slide 15; `Lecture05`, Slide 68).
   - I residui d'interfaccia responsabili del legame tra Myotrophin e Actin-Capping Protein (identificati oggettivamente per sovrapposizione strutturale con il complesso cristallografico `3AAA` entro un cutoff di contatto atomico $\le 4.5\text{ Å}$) devono essere documentati e **mascherati (esclusi)** da qualsiasi mutazione.

### B. Scelta del Force Field e Setup della Dinamica Molecolare
1. **Giustificazione del Force Field (`CHARMM36m`):**
   - Utilizzare **CHARMM36m** (`top_all36_prot.rtf` + `par_all36m_prot.prm` + `toppar_water_ions.str` da `Tutorials/charmm36/`), lo standard stato-dell'arte per proteine globulari e a repeat, che corregge la sovra-stabilizzazione elicoidale delle versioni precedenti (`CHARMM22/27`).
2. **Solvatazione e Forza Ionica Fisiologica:**
   - Per simulazioni in solvente esplicito (`TIP3P`), garantire un padding minimo che impedisca alla proteina di "vedere" la propria immagine periodica attraverso le Periodic Boundary Conditions (**cutoff non-bonded $12.0\text{ Å}$** $\rightarrow$ distanza minima tra immagini $> 24\text{ Å}$, ovvero padding $\ge 10\text{--}12\text{ Å}$ per lato; `Lecture03`, Slide 36–39).
   - Neutralizzare la carica netta del sistema e impostare una concentrazione salina fisiologica (**$0.15\text{ M NaCl}$**) tramite `autoionize` (`Lecture03`, Slide 39).
3. **Protocollo "Matched Comparison" e Repliche Multi-Seed:**
   - Il confronto tra Wild-Type (`WT`), mutanti singoli (`Mut1`, `Mut2`) e doppio mutante (`Double`) deve avvenire a **condizioni rigorosamente identiche** (*Matched Protocol*, `Laboratory_Project_2026.pdf`, Slide 15 & Slide 21).
   - Ogni sistema deve essere simulato con **almeno $n = 3$ repliche indipendenti** generate variando il `seed` del generatore di numeri casuali per l'assegnazione delle velocità iniziali di Maxwell-Boltzmann (`Lecture04`, Slide 28–30; *Bekker et al., 2019*).

### C. Analisi delle Traiettorie e Rigore Statistico
1. **Pre-processing Obbligatorio della Traiettoria (`pbc wrap` + `measure fit`):**
   - Prima di calcolare qualsiasi metrica posizionale (`RMSD`, `RMSF`), è obbligatorio ricomporre le condizioni periodiche al contorno (`pbc wrap -compound residue -center com -centersel "protein" -all`) e rimuovere i moti roto-traslazionali globali allineando tutti i frame sulla struttura iniziale (`frame 0`) tramite `measure fit` (`Laboratory_Project_2026.pdf`, Slide 22; `Lecture04`, Slide 31–34).
2. **Statistica e Incertezza:**
   - Riportare sempre **Media $\pm$ Deviazione Standard ($\mu \pm \sigma$)** calcolata tra le repliche indipendenti sulla finestra temporale equilibrata finale, verificando la significatività statistica delle differenze osservate.

### D. Compatibilità Multipiattaforma (Linux + Windows), Visualizzazioni, Cartella `STATUS/` e Gestione Git
1. **📊 Regola Formato Grafici (`.pdf` Vettoriale Obbligatorio):**
   - Tutti i grafici 2D generati in Python (`matplotlib`/`seaborn`) devono essere salvati obbligatoriamente in **formato vettoriale `.pdf`** (`format="pdf"`, non `.png`) dentro `06_presentation/`, per garantire qualità tipografica perfetta a qualsiasi livello di zoom nelle slide e nei report.
2. **📓 Tracciamento Continuo nella Cartella `STATUS/` (`discovery.md` e `SESSION<N>.md`):**
   - **`STATUS/discovery.md`**: Registro cumulativo di tutte le scoperte scientifiche, biofisiche e HPC emerse durante il lavoro (es. reti di H-bond delle istidine, deviazioni dal consenso Ankyrin, residui funzionali, speedup AVX-512 su core fisici vs SMT, ecc.). Va aggiornato ogni volta che emerge un risultato degno di nota per la presentazione.
   - **`STATUS/SESSION0.md`, `SESSION1.md`, ...**: Diario dettagliato per ogni sessione di lavoro che documenta cosa è stato fatto, gli script creati/validati, i risultati ottenuti e la scaletta operativa per la sessione successiva.
3. **📸 Regola "Momento Foto per le Slide" (Immagini VMD ad Alta Risoluzione):**
   - In ogni fase del progetto (`ROAD1.md`, `ROAD2.md`, ...) in cui si ottiene un risultato strutturale rilevante per la presentazione finale, deve essere sempre presente un box dedicato **`📸 MOMENTO FOTO PER LE SLIDE (VMD)`** (accompagnato da uno script `.tcl` di visualizzazione in `scripts/viz_*.tcl`).
   - Lo script imposta automaticamente in VMD le buone norme grafiche da pubblicazione/presentazione (`color Display Background white`, `axes location Off`, `display projection Orthographic`, `display depthcue off`, `NewCartoon` + `Licorice`) e indica all'utente come salvare l'immagine in `06_presentation/` tramite **`File -> Render -> Snapshot (o Tachyon)`**.
4. **Percorsi Relativi:** Tutti gli script `.tcl`, `.conf` e `.py` devono utilizzare esclusivamente **percorsi relativi** con slash `/` (supportati nativamente da Tcl/VMD, NAMD e Python sia su Linux che su Windows).
5. **Igiene del Repository Git:** Non committare mai file binari di traiettoria o checkpoint pesanti (`.dcd`, `.coor`, `.vel`, `.xsc`, `.bak`, `.old`). Committare sempre script, file di configurazione, strutture `.pdb`/`.psf` di riferimento, tabelle di output `.dat`/`.csv`/`.fxout`, figure `.pdf` in `06_presentation/`, diario in `STATUS/` e documentazione in `ROADMAP/`.
