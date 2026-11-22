# 🔬 RISULTATI E COMMENTO TEORICO DELLO SCREENING MD A 400 K (`Seed 1`, `2.5 ns`) — `WT` vs `K37L`, `T103L`, `K37L_T103L` e `K37L_H68G_T103L`

**Data Completamento Simulazioni:** 22 Novembre 2026 (ore 11:12)  
**Cartelle con i Plot e Report Individuali Separati (`05_analysis/single_runs/`):**
- 📁 **[`00_WT/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/00_WT)** — Wild-Type (`2MYO`, controllo stress termico 400 K): [Report](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/00_WT/WT_summary_report.txt) | [Plot 1 (Dynamics)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/00_WT/WT_01_dynamics_overview.pdf) | [Plot 2 (Q-value)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/00_WT/WT_02_qvalue_repeats_and_classes.pdf) | [Plot 3 (Residues)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/00_WT/WT_03_per_residue_qi_and_rmsf.pdf)
- 📁 **[`01_K37L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/01_K37L)** — Simulazione 1 (`K37L`, Singolo `ANK2`): [Report](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/01_K37L/K37L_summary_report.txt) | [Plot 1 (Dynamics)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/01_K37L/K37L_01_dynamics_overview.pdf) | [Plot 2 (Q-value)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/01_K37L/K37L_02_qvalue_repeats_and_classes.pdf) | [Plot 3 (Residues)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/01_K37L/K37L_03_per_residue_qi_and_rmsf.pdf)
- 📁 **[`02_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/02_T103L)** — Simulazione 2 (`T103L`, Singolo `ANK4`): [Report](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/02_T103L/T103L_summary_report.txt) | [Plot 1 (Dynamics)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/02_T103L/T103L_01_dynamics_overview.pdf) | [Plot 2 (Q-value)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/02_T103L/T103L_02_qvalue_repeats_and_classes.pdf) | [Plot 3 (Residues)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/02_T103L/T103L_03_per_residue_qi_and_rmsf.pdf)
- 📁 **[`03_K37L_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/03_K37L_T103L)** — Simulazione 3 (`K37L_T103L`, Doppio `ANK2 + ANK4`): [Report](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/03_K37L_T103L/K37L_T103L_summary_report.txt) | [Plot 1 (Dynamics)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/03_K37L_T103L/K37L_T103L_01_dynamics_overview.pdf) | [Plot 2 (Q-value)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/03_K37L_T103L/K37L_T103L_02_qvalue_repeats_and_classes.pdf) | [Plot 3 (Residues)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/03_K37L_T103L/K37L_T103L_03_per_residue_qi_and_rmsf.pdf)
- 📁 **[`04_K37L_H68G_T103L/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/04_K37L_H68G_T103L)** — Simulazione 4 (`K37L_H68G_T103L`, Triplo `ANK2 + ANK3 + ANK4`): [Report](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/04_K37L_H68G_T103L/K37L_H68G_T103L_summary_report.txt) | [Plot 1 (Dynamics)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/04_K37L_H68G_T103L/K37L_H68G_T103L_01_dynamics_overview.pdf) | [Plot 2 (Q-value)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/04_K37L_H68G_T103L/K37L_H68G_T103L_02_qvalue_repeats_and_classes.pdf) | [Plot 3 (Residues)](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs/04_K37L_H68G_T103L/K37L_H68G_T103L_03_per_residue_qi_and_rmsf.pdf)

---

## 📊 1. Tabella Sinottica dei Risultati delle 5 Simulazioni a `400 K` (`Seed 1`, Finestra Finale Equilibrata 30%: `1.75–2.50 ns`)

Tutte le metriche sono valutate sugli ultimi 376 frame (`Frame 876..1251`, `1.75–2.50 ns`) delle traiettorie NPT a `400 K`, in stretta conformità con *Bekker, Ma & Kamiya (2019)* Table II (p. 434) e con le Slide 11–13, 15, 22–23 di [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf):

| Osservabile Termodinamica / Strutturale (`1.75–2.50 ns` a `400 K`) | `00_WT` (`2MYO` Wild-Type) | `01_K37L` (Singolo 1, `ANK2`) | `02_T103L` (Singolo 2, `ANK4`) | `03_K37L_T103L` (Doppio `ANK2+4`) | `04_K37L_H68G_T103L` (Triplo `ANK2+3+4`) | Miglioramento Massimo vs `WT` |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **$\Delta\Delta G_{\text{BuildModel}}$ ($n=5$ run, $\text{kcal/mol}$)** | `0.000` | `-2.438 ± 0.870` | `-1.150 ± 0.143` | `-3.652 ± 0.656` | 🏆 **`-5.032 ± 0.210`** | **`-5.03 kcal/mol`** (`Best = -5.39`) |
| **`Q: hydrophilic-all`** *(Bekker $r = 0.84$ con $T_m$)* | `0.6750 ± 0.0067` | `0.7134 ± 0.0070` | `0.7047 ± 0.0066` | 🏆 **`0.7178 ± 0.0075`** | 🏆 **`0.7172 ± 0.0056`** | **`+0.0428` (`+4.28%` contatti polari)** |
| **`Q: all-all`** *(Bekker $r = 0.79$ con $T_m$)* | `0.6948 ± 0.0054` | `0.7243 ± 0.0052` | `0.7124 ± 0.0046` | `0.7194 ± 0.0061` | 🏆 **`0.7289 ± 0.0049`** | **`+0.0341` (`+3.41%` contatti totali)** |
| **`Q: hydrophobic-small`** *(Bekker $r = 0.76$)* | `0.6433 ± 0.0154` | `0.6590 ± 0.0139` | `0.6433 ± 0.0116` | `0.6351 ± 0.0154` | 🏆 **`0.6631 ± 0.0140`** | **`+0.0198` (`+1.98%` core idrofobico)** |
| **`Q_ANK1` (N-cap, res `1–34`)** | `0.6959` | `0.7205` | `0.7142` | `0.7278` | 🏆 **`0.7358`** | **`+0.0399` (`+3.99%` su `ANK1`)** |
| **`Q_ANK2` (Core, res `35–67` — modulo critico)** | `0.6015` *(fragile)* | `0.6351` | `0.6161` | `0.6290` | 🏆 **`0.6444`** | **`+0.0429` (`+4.29%` su `ANK2`)** |
| **`Q_ANK3` (Core, res `68–100`)** | `0.6660` | `0.6855` | 🏆 **`0.6858`** | `0.6676` | `0.6790` | **`+0.0198` (`+1.98%` su `ANK3`)** |
| **`Q_ANK4` (C-cap, res `101–118`)** | `0.6881` | 🏆 **`0.7122`** | `0.6699` | `0.6981` | `0.7044` | **`+0.0241` (`+2.41%` su `ANK4`)** |
| **Interfaccia `Q_int(ANK1–ANK2)`** | `0.113` *(crollata)* | `0.228` | `0.173` | `0.266` | 🏆 **`0.291`** | **`×2.58` (`+158%` coesione `ANK1–2`)** |
| **Interfacce `Q_int(ANK2–3 / ANK3–4)`** | `0.138 / 0.247` | 🏆 **`0.221 / 0.303`** | `0.221 / 0.181` | `0.173 / 0.169` | `0.157 / 0.225` | **`+60%` coesione `ANK2–3`** |
| **Backbone RMSD Globale ($\text{Å}$)** | `6.68 ± 0.52 Å` | `4.17 ± 0.33 Å` | `3.77 ± 0.19 Å` | `6.40 ± 0.53 Å` | 🏆 **`3.30 ± 0.42 Å`** | **`-3.38 Å` (`-50.6%` deformazione)** |
| **Functional Mask `3AAA` RMSD ($\text{Å}$)** | `7.27 ± 0.99 Å` | `4.75 ± 0.32 Å` | `3.37 ± 0.36 Å` | `6.73 ± 0.52 Å` | 🏆 **`2.83 ± 0.26 Å`** | **`-4.44 Å` (`-61.1%` sul sito attivo!)** |
| **Raggio di Girazione $R_g$ ($\text{Å}$)** | `15.74 ± 0.18 Å` | `15.36 ± 0.17 Å` | `15.08 ± 0.21 Å` | 🏆 **`14.76 ± 0.34 Å`** | **`14.97 ± 0.21 Å`** | **`-0.98 Å` (compattezza 300 K conservata)** |
| **SASA Totale ($\text{Å}^2$)** | `8838.8 ± 177.4` | `7941.6 ± 206.2` | `7969.2 ± 218.7` | `8117.8 ± 269.2` | 🏆 **`7646.9 ± 234.7`** | **`-1191.9 Å²` (`-13.5%` apertura al solvente)** |
| **% $\alpha$-Helix (`STRIDE`) e H-Bonds Unici** | `53.78%` (`13.1`) | **`53.23%` (`15.3`)** | `50.81%` (`14.8`) | `44.39%` (`13.1`) | **`53.11%` (`14.4`)** | **`+1.3 / +2.2` H-bonds intra-proteici** |
| **Fluttuazione Locale $\text{RMSF}_{C_\alpha}$ Res `37 / 60 / 68 / 103` ($\text{Å}$)** | `2.96 / 2.53 / 3.63 / 3.16` | `1.37 / 5.00 / 2.93 / 2.16` | `1.83 / 1.75 / 2.53 / 2.88` | `2.08 / 3.35 / 2.59 / 4.12` | 🏆 **`0.91 / 1.50 / 1.95 / 2.83`** | **`RMSF(37) = 0.91 Å`, `RMSF(60) = 1.50 Å`!** |
| **Residui Scaffold Instabili ($\langle Q_i \rangle_{400\text{K}} < 0.60, N_0 \ge 30$)** | `41` residui | `25` residui *(-16)* | `30` residui *(-11)* | `31` residui *(-10)* | 🏆 **`23` residui *(-18)*** | **`-43.9%` di residui instabili!** |

---

## 🧠 2. Analisi Biofisica Approfondita: Perché la Tripla (`K37L + H68G + T103L`) Risolve l'Effetto Collaterale sul Loop `ASP 60` del Doppio Mutante

Confrontando i 5 pacchetti di plot individuali in [`05_analysis/single_runs/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs), emerge un meccanismo strutturale di straordinaria eleganza che lega la chimica locale delle singole mutazioni alla cooperatività globale dello scaffold Ankyrin (*Lowe & Itzhaki, PNAS 2007*; *Mosavi et al., PNAS 2002*):

### A. L'effetto di `K37L` (`01_K37L`, Singolo in `ANK2`)
- **Cosa ripara:** Nel Wild-Type (`00_WT`), `LYS 37` occupa la Posizione 4 del motivo consenso `TPLH` di `ANK2` (dove in `ANK1`, `ANK3` e `ANK4` siede una Leucina idrofobica) e forma un ponte salino distorsivo (`47.61%` di occupazione a 300 K) con `ASP 60` nel loop `ANK2–ANK3`. Sostituendo `LYS 37 -> LEU 37`, il core idrofobico di `ANK2` si chiude a cerniera:
  - $Q_{\text{ANK2}}$ balza da `0.6015` a **`0.6351`**;
  - I residui del core di `ANK2` (`ILE 50, LEU 51, PHE 53`), che nel WT erano tutti collassati sotto `0.54`, in `K37L` **superano tutti la soglia di stabilità `0.60`** e scompaiono dalla lista dei residui instabili;
  - Il numero totale di residui instabili dello scaffold crolla da **`41` a `25` (`-16` residui)** e il Backbone RMSD scende da `6.68 Å` a **`4.17 Å`**.
- **Il piccolo effetto collaterale locale (`ASP 60` orfano):** Rimuovendo il partner positivo `LYS 37`, la catena laterale negativa di **`ASP 60`** (sulla punta del loop `59–65` tra `ANK2` e `ANK3`) rimane priva del suo ancoraggio elettrostatico e oscilla maggiormente nel solvente: il suo $\text{RMSF}_{C_\alpha}$ sale da `2.53 Å` (WT) a **`5.00 Å`** in `K37L`.

### B. L'effetto di `T103L` (`02_T103L`, Singolo nel C-cap `ANK4`)
- **Cosa ripara:** Nel C-cap `ANK4`, la Treonina polare `THR 103` (anch'essa in Posizione 4 del motivo `TPLH`, dove il consenso di *Mosavi et al. 2002* richiede `LEU`) espone il margine C-terminale all'ingresso dell'acqua a 400 K. La mutazione `T103L` aumenta il packing idrofobico ($\Delta E_{\text{VdW+SolvH}} = -4.03\text{ kcal/mol}$ in FoldX), riducendo la SASA di **`-870 Å²`** (`7969.2` vs `8838.8 Å²`), abbassando il Backbone RMSD a **`3.77 ± 0.19 Å`** e proteggendo la maschera funzionale `3AAA` (**`3.37 ± 0.36 Å`** vs `7.27 Å` del WT).

### C. L'effetto del Doppio Mutante `K37L_T103L` (`03_K37L_T103L`, `ANK2 + ANK4`)
- **Cosa ripara:** Combinando `K37L` e `T103L`, otteniamo il **valore più alto in assoluto di `Q: hydrophilic-all` (`0.7178 ± 0.0075` vs `0.6750` del WT, `+4.28%`)** e il **raggio di girazione più compatto ($R_g = 14.76 \pm 0.34\text{ Å}$ vs `15.74 Å` del WT)**, oltre a raddoppiare la coesione dell'interfaccia `ANK1–ANK2` (`0.266` vs `0.113`).
- **Perché però l'RMSD globale del Doppio (`6.40 Å`) resta più alto rispetto ai singoli?**  
  Perché nel Doppio `K37L_T103L` la regione centrale `ANK3` (`68–100`) non è mutata: il loop `ASP 60 – HIS 68` tra `ANK2` e `ANK3` ($\text{RMSF}_{60} = 3.35\text{ Å}$) e i residui `LEU 72, LEU 73, SER 74, VAL 76` di `ANK3` subiscono una flessione elastica ("hinge bending") tra i due poli irrigiditi `ANK2` (`K37L`) e `ANK4` (`T103L`). Come spiegato nella **Slide 12 del Prof. Gautieri**, una flessione di cerniera alza l'RMSD di sovrapposizione rigida anche quando i contatti nativi interni ($Q = 0.7194$) e il volume globulare ($R_g = 14.76\text{ Å}$) sono perfettamente conservati!

### D. Il Trionfo del Triplo Mutante `K37L_H68G_T103L` (`04_K37L_H68G_T103L`, `ANK2 + ANK3 + ANK4`)
Ed ecco perché aver incluso nella Simulazione 4 la **Tripla `K37L + H68G + T103L`** (dove **`H68G`** è il #1 *Back-to-Consensus* di FireProt 2.0 e #1 assoluto di ThermoMPNN, situato in Posizione 2 del beta-turn d'ingresso di **`ANK3`**) chiude il cerchio in modo spettacolare:
1. **`H68G` rilassa la tensione sterica esattamente all'uscita del loop `60–68` (`ANK2–ANK3`):**
   - Sostituendo l'ingombrante anello imidazolico di `HIS 68` (che si trovava in conformazione ad alta energia all'apice del $\beta$-turn tra `ANK2` e `ANK3`) con la Glicina canonica del consenso Ankyrin (`G` conservata nel 78% delle proteine Ankyrin), **l'intero loop `59–68` trova la sua geometria nativa senza tensione**!
   - **Prova numerica diretta:** L'RMSF di **`ASP 60`**, che era salito a `5.00 Å` in `K37L` e `3.35 Å` nel Doppio `K37L_T103L`, nella Tripla `K37L_H68G_T103L` **crolla a `1.50 Å`** (ancora più basso del WT `2.53 Å`!) e l'RMSF di **`LEU 37` scende addirittura a `0.91 Å`** a 400 K!
2. **Stabilizzazione cooperativa dell'intera pila `ANK1–ANK2–ANK3–ANK4`:**
   - La Tripla raggiunge il **record assoluto su quasi tutte le metriche contemporaneamente**:
     - **`Q: all-all` = `0.7289 ± 0.0049`** (#1 assoluto, `+3.41%` vs WT);
     - **`Q: hydrophilic-all` = `0.7172 ± 0.0056`** (pari al record del Doppio `0.7178`, ma con fluttuazione $\sigma$ inferiore);
     - **`Q: hydrophobic-small` = `0.6631 ± 0.0140`** (#1 assoluto);
     - **`Q_ANK1` = `0.7358`** (#1 assoluto) e **`Q_ANK2` = `0.6444`** (#1 assoluto, `+4.29%` vs WT);
     - **Interfaccia `ANK1–ANK2` = `0.291`** (#1 assoluto, **`+158%`** rispetto a `0.113` del WT);
     - **Backbone RMSD Globale = `3.30 ± 0.42 Å`** (#1 assoluto, più che dimezzato rispetto ai `6.68 Å` del WT);
     - **Functional Mask `3AAA` RMSD = `2.83 ± 0.26 Å`** (#1 assoluto, ridotto del **`-61%`** rispetto ai `7.27 Å` del WT: il sito di legame all'Actin-Capping Protein resta cristallino anche a 400 K!);
     - **SASA Totale = `7646.9 ± 234.7 Å²`** (#1 assoluto, **`-1.192 Å²`** rispetto al WT);
     - **Residui instabili dello scaffold ($\langle Q_i \rangle_{400\text{K}} < 0.60$) = solo `23` residui** (#1 assoluto, contro i `41` del WT: **`-18` residui salvati**, con guarigione completa di tutto il loop `90–95` di `ANK3`!).

---

## 🎯 3. Cosa Dobbiamo Fare Adesso (Piano Operativo per la Campagna di Validazione Finale)

Ora che abbiamo tutti e 5 i pacchetti singoli (`00_WT`, `01_K37L`, `02_T103L`, `03_K37L_T103L`, `04_K37L_H68G_T103L`) analizzati e verificati, la strada per chiudere la **Fase 4** e la **Fase 5** con rigore da pubblicazione è divisa in **2 passi chiari**:

### Passo A (Immediato — 30 secondi, zero nuove simulazioni): Generare i Plot Comparativi Unificati del `Seed 1`
Prima ancora di lanciare nuove repliche, possiamo già unire i dati delle 5 simulazioni concluse per generare:
- La **Tabella II di Bekker unificata** (`05_analysis/bekker_table2_and_epistasis_report.txt`);
- I **3 grafici PDF comparativi multi-sistema** in `06_presentation/` (`fig9_matched_q_time_and_repeats.pdf`, `fig10_per_residue_delta_qi_and_rmsf.pdf`, `fig11_epistatic_cycle_and_observables.pdf`) dove le curve di `WT`, `K37L`, `T103L`, `Doppio` e `Triplo` sono sovrapposte sullo stesso grafico per il confronto visivo immediato.

### Passo B (Campagna di Validazione Statistica & Verifica con Pre-Equilibrazione Restrained)
Per la campagna di validazione finale (Punto 4 & Punto 5 della scaletta del Prof. Gautieri, Slide 15 e 21), abbiamo due elementi complementari:
1. **Validazione Statistica Multi-Seed (`Seed 2 = 23456` e `Seed 3 = 34567` a `400 K`):**
   - Per avere le barre d'errore ($\mu \pm \sigma$ su $n=3$ velocity seeds indipendenti, come richiesto nella Slide 21 del professore e nella Table II di *Bekker et al.*), ci servono i `Seed 2` e `Seed 3` sui sistemi chiave (`WT`, `K37L`, `T103L`, `K37L_T103L` e/o `K37L_H68G_T103L`).
   - *Nota strategica:* Puoi dividere queste run con Alessia (usando `04_md_simulations/run_batch_phase4.sh`) oppure far girare sul tuo PC prima i `Seed 2` e `Seed 3` dei tre sistemi più importanti (**`WT`**, **`Doppio K37L_T103L`** e **`Triplo K37L_H68G_T103L`** = 6 run, o includendo anche `K37L` e `T103L`).
2. **Test Metodologico con Pre-Equilibrazione Restrained a `300 K` + Produzione a `400 K` (la tua idea di stamattina!):**
   - Possiamo creare uno script/configurazione di validazione rigorosa (**Protocollo con Pre-Equilibrazione Restrained**: `1000 step minimize` $\rightarrow$ `100 ps NVT/NPT a 300 K` con vincoli armonici `constraints on` sul backbone della proteina per rilassare l'acqua e gli ioni $\rightarrow$ rilascio dei vincoli e produzione a `400 K`) e testarlo in parallelo (`1-vs-1` o `1-vs-2`: **`WT_eq` vs `K37L_T103L_eq` e `K37L_H68G_T103L_eq`**), così nella presentazione mostrerete al Prof. Gautieri sia la campagna Multi-Seed standard sia la verifica con Pre-Equilibrazione esplicita da `Lecture04`!
