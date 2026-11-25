# 🎤 DISCORSO_PROF_RISULTATI_INTERMEDI.md — Scaletta e Discorso Pronto per Presentare i Risultati Intermedi al Prof. Gautieri

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Proteina Assegnata:** Myotrophin (*Rattus norvegicus*, PDB ID: [`2MYO`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/2MYO.pdb), 118 aa — 4 Ankyrin Repeats)

---

## 🛠️ PRIMA DI INIZIARE: I 3 Comandi da Rilanciare Passo-Passo (per rigenerare tutti i $Q$ che partono da `1.0` e stanno a `0.80–0.85`)

Dalla cartella `/home/olly/UNI/BioMolEng2026/laboratorio/BML`, lancia questi **3 comandi in ordine**:

### Passo 1 — Ricalcolo Baseline Wild-Type (`300 K` vs `400 K`) e Plot di Fase 2 (inclusi i plot di Alessia in `report/02_utilities/`):
```bash
vmd -dispdev text -e scripts/03_analyze_wt_baseline.tcl && \
python3 scripts/plot_wt_baseline.py && \
python3 report/02_scripts/02_plot_rmsd_q_wt.py && \
python3 report/02_scripts/03_plot_wt_remaining.py && \
python3 report/02_scripts/04_plot_wt_structural_regions.py
```

### Passo 2 — Ricalcolo delle 5 Simulazioni Comparative a `400 K` (`WT`, `K37L`, `T103L`, `Doppio`, `Triplo`) e Plot Singoli Separati:
```bash
vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl && \
python3 scripts/10_plot_single_runs_400K.py
```

### Passo 3 — Generazione della Tabella II di Bekker Unificata e dei 3 Grafici Comparativi (`fig9`, `fig10`, `fig11`):
```bash
python3 scripts/10_plot_final_validation_bekker.py && \
cat 05_analysis/bekker_table2_and_epistasis_report.txt
```

---

## 🎙️ DISCORSO PRONTO PER IL PROF. GAUTIERI (Diviso nei 5 Punti della sua Scaletta — Slide 15)

Ecco la traccia discorsiva esatta, punto per punto, con i file/grafici da mostrargli sullo schermo mentre parli:

---

### 🔹 1. Analisi Strutturale e Funzionale di `2MYO` (Punto 1 della Scaletta)
> **Cosa mostrare a schermo:** [`report/01_protein_characterization.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/01_protein_characterization.md) oppure [`06_presentation/fig6_ankyrin_consensus_alignment.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig6_ankyrin_consensus_alignment.pdf)
>
> **Cosa dire al Professore:**
> *"Professore, la nostra proteina è la **Myotrophin (`2MYO`)**, una proteina intracellulare di 118 amminoacidi formata da **4 Ankyrin Repeats** (`ANK1` N-cap, due repeat centrali completi `ANK2` e `ANK3`, e `ANK4` C-cap), la cui funzione fisiologica principale è legare e inibire stericamente l'**Actin-Capping Protein (CP)** per regolare la polimerizzazione dell'actina.*
> *Abbiamo verificato che la struttura NMR `2MYO` è completa (118 residui continui, zero atomi pesanti mancanti e 3 cisteine tutte ridotte con distanza $S_\gamma\text{--}S_\gamma > 9.4\text{ Å}$). Prima di costruire la topologia `CHARMM36m`, abbiamo assegnato lo stato di protonazione delle 4 istidine ispezionando la rete di legami idrogeno locali (`HIS 40, 68, 80 -> HSD`, mentre `HIS 67 -> HSE` per evitare repulsione con `LYS 37` e donare un H-bond a `GLY 35 O`).*
> *Soprattutto, per non alterare la funzione biologica durante l'ingegnerizzazione, abbiamo sovrapposto `2MYO` al complesso cristallografico Myotrophin–Actin Capping Protein (**PDB `3AAA`**, $\text{RMSD} = 2.36\text{ Å}$) e abbiamo bloccato in una **Functional Mask intoccabile tutti i 24 residui a distanza $< 4.5\text{ Å}$ dal partner di legame**."*

---

### 🔹 2. Caratterizzazione Baseline del Wild-Type a `300 K` e `400 K` (Punto 2 della Scaletta)
> **Cosa mostrare a schermo:** [`report/02_utilities/02_WT_RMSD_time.png`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_utilities/02_WT_RMSD_time.png), [`report/02_utilities/03_WT_Q_time.png`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_utilities/03_WT_Q_time.png) e [`06_presentation/fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf)
>
> **Cosa dire al Professore:**
> *"Abbiamo solvatato il sistema in acqua esplicita `TIP3P` (padding $12\text{ Å}$) e neutralizzato a $0.15\text{ M NaCl}$ (`21.461` atomi, carica netta zero). Dopo aver fatto il benchmark su NAMD 3 (~`12.7 ns/day` su 8 core fisici), abbiamo eseguito le due simulazioni pilota del Wild-Type: `1.5 ns` a **`300 K`** come controllo fisiologico e `2.5 ns` a **`400 K`** sotto stress termico, seguendo l'approccio di Bekker et al. (2019).*
> *Calcolando la frazione di contatti nativi $Q(t)$ rispetto al primo frame della traiettoria minimizzata (`Frame 0` del `.dcd`, così da partire esattamente da $Q(0) = 1.0$ senza artefatti da minimizzazione del modello NMR):*
> - *A **`300 K`** la proteina è molto stabile: l'RMSD si stabilizza a `~2.6 Å`, il raggio di girazione resta `14.5 Å` e il **$Q$-value globale rimane alto intorno a `0.85–0.88`**.*
> - *A **`400 K`** osserviamo un'espansione strutturale ($R_g$ sale a `15.7 Å`, SASA sale a `8840 Å²`, RMSD supera i `6.6 Å`) e il **$Q$-value globale scende nella fascia `0.80–0.84`**.*
> *Analizzando l'RMSF e il $Q_i$ per singolo residuo negli ultimi 30% della traiettoria a 400 K, abbiamo identificato **3 punti deboli strutturali ben precisi** fuori dalla maschera funzionale:*
> 1. *Il **primo repeat centrale (`ANK2`, residui 35–67)** è il modulo che perde più contatti nativi, innescato da **`LYS 37`** (che mostra un picco di RMSF da `0.95` a `2.96 Å` e forma un ponte salino distorsivo con `ASP 60`, destabilizzando il core idrofobico vicino `ILE 50 – LEU 51`);*
> 2. *La **giunzione centrale tra i due repeat (`residui 60–69`, attorno a `ASP 60 – HIS 68`)**, dove `HIS 68` è anche un outlier di Ramachandran nel `REMARK 500`;*
> 3. *L'**ingresso del C-cap (`ANK4`, `LEU 102 – THR 103`)**, dove si perde il legame a idrogeno di capping `THR103 – ASP100` (che scende dal `15.6%` a 300 K allo `0.3%` a 400 K)."*

---

### 🔹 3. Ipotesi Mutazionali, Screening Multi-Tool e `FoldX BuildModel` (Punto 3 della Scaletta)
> **Cosa mostrare a schermo:** [`06_presentation/fig8_consensus_energy_breakdown.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8_consensus_energy_breakdown.pdf) e [`06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf)
>
> **Cosa dire al Professore:**
> *"Per scegliere le mutazioni non ci siamo fidati di un singolo tool, ma abbiamo incrociato i punti deboli dinamici della MD a 400 K con **4 predittori indipendenti** sulle 88 posizioni non mascherate:*
> 1. ***FoldX 5** (`RepairPDB` $\rightarrow$ `PositionScan` su tutte le $88 \times 20 = 1.760$ mutazioni $\rightarrow$ `BuildModel` con 5 run e riferimenti WT accoppiati);*
> 2. ***FireProt 2.0** (ibrido energetico Rosetta + evolutivo Back-to-Consensus);*
> 3. ***ThermoMPNN** (Graph Neural Network addestrata su Megascale $T_m$);*
> 4. *L'**allineamento al consenso universale degli Ankyrin Repeats** di Mosavi et al. (2002).*
>
> *Da questo filtro sono emerse **3 mutazioni chiave**, tutte e tre di tipo **Back-to-Consensus** e ciascuna mirata a uno dei 3 punti deboli visti nella MD:*
> - ***`K37L` (in `ANK2`):** Sostituisce la Lisina 37 con la **Leucina canonica del motivo `TPLH`** (già presente nella stessa posizione in `ANK1`, `ANK3` e `ANK4`). In FoldX `BuildModel` dà **$\Delta\Delta G = -2.44 \pm 0.87\text{ kcal/mol}$** ed è confermata sia da FireProt che da ThermoMPNN.*
> - ***`T103L` (nel C-cap `ANK4`):** Ripristina anche qui la **Leucina canonica in posizione 4 di `ANK4`**, migliorando l'impacchettamento idrofobico ($\Delta E_{\text{VdW+SolvH}} = -4.03\text{ kcal/mol}$, $\Delta\Delta G_{\text{BuildModel}} = -1.15 \pm 0.14\text{ kcal/mol}$).*
> - ***`H68G` (alla giunzione `ANK2–ANK3`):** È la **#1 scelta Back-to-Consensus di FireProt 2.0** (conservata nel 78% delle proteine Ankyrin) e la **#1 assoluta di ThermoMPNN** (`-1.44 kcal/mol`), con $\Delta\Delta G_{\text{BuildModel}} = -1.73 \pm 0.01\text{ kcal/mol}$: sostituisce l'istidina 68 in tensione sterica sul $\beta$-turn con una Glicina.*
> - *Quando le abbiamo combinate in `FoldX BuildModel`, il **Doppio Mutante `K37L + T103L`** è risultato perfettamente additivo (**$-3.65 \pm 0.66\text{ kcal/mol}$**) e il **Triplo Mutante `K37L + H68G + T103L`** ha raggiunto **$-5.03 \pm 0.21\text{ kcal/mol}$**."*

---

### 🔹 4. Risultati Intermedi dello Screening MD a `400 K` (`WT` vs Singoli, Doppio e Triplo) (Punto 4 & 5)
> **Cosa mostrare a schermo:** [`06_presentation/fig9_matched_q_time_and_repeats.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig9_matched_q_time_and_repeats.pdf), [`06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf) e [`06_presentation/fig11_epistatic_cycle_and_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig11_epistatic_cycle_and_observables.pdf)
>
> **Cosa dire al Professore:**
> *"Abbiamo quindi costruito con lo stesso identico protocollo (`CHARMM36m + TIP3P 12 Å + 0.15 M NaCl`, stessa box PBC e carica zero) i sistemi solvatati per i due singoli (`K37L`, `T103L`), il Doppio (`K37L+T103L`) e il Triplo (`K37L+H68G+T103L`) e abbiamo eseguito un primo ciclo di screening dinamico a **`400 K` per `2.5 ns`** a parità di velocity seed.*
> *I risultati dinamici confermano il design e rivelano un effetto cooperativo molto interessante:*
> 1. ***Tutti i mutanti migliorano il $Q$-value (`all-all` e `hydrophilic-all`), la compattezza ($R_g$) e riducono la SASA rispetto al Wild-Type.***
> 2. *In particolare, **`K37L`** stabilizza il core idrofobico di `ANK2` (`ILE 50, LEU 51, PHE 53` recuperano i loro contatti nativi) e dimezza l'RMSD globale da `6.68 Å` a `4.17 Å`, ma rompendo il ponte salino `LYS 37 – ASP 60` lascia il loop centrale `ASP 60` più mobile ($\text{RMSF}_{60} = 5.00\text{ Å}$ in `K37L` e `3.35 Å` nel Doppio `K37L+T103L`).*
> 3. *Quando aggiungiamo **`H68G` nel Triplo Mutante (`K37L + H68G + T103L`)**, la Glicina 68 rilassa la cerniera tra `ANK2` e `ANK3` e **abbatte la fluttuazione del loop `ASP 60` a soli `1.50 Å`** (e quella di `LEU 37` a **`0.91 Å`**!). Di conseguenza, il **Triplo Mutante** ottiene il miglior **$Q$-value globale**, dimezza il **Backbone RMSD globale a `3.30 Å`** (contro `6.68 Å` del WT), protegge il **sito di legame all'actina (`3AAA` RMSD = `2.83 Å` contro `7.27 Å` del WT)** e riduce del **44%** il numero di residui instabili!"*

---

### 🔹 5. Prossimo Passo: Validazione Finale Multi-Seed con Pre-Equilibrazione Esplicita a `300 K`
> **Cosa dire al Professore per chiudere:**
> *"Come ultimo passo di validazione rigorosa prima della presentazione finale, avendo notato che nelle run pilota la dinamica NPT partiva subito dopo la minimizzazione, stiamo ora lanciando la **campagna di validazione multi-seed ($n=3$ velocity seeds indipendenti)** sul confronto diretto tra **Wild-Type** e i nostri candidati migliori (**Doppio/Triplo Mutante**) inserendo anche uno stadio esplicito di **Pre-Equilibrazione a `300 K` con position restraints sul backbone (`100 ps`)** prima dello stress termico a `400 K`, così da presentare le medie $\mu \pm \sigma$ definitive."*
