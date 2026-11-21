# 🎯 MUTATION_SELECTION.md — Criteri Motivati, Processo Logico e Risultati Finali della Scelta delle Mutazioni

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Target:** Myotrophin (`2MYO`, 118 amminoacidi)  
**Regola Fondamentale ([`AGENTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/AGENTS.md)):** *Zero scelte arbitrarie. Ogni inclusione o esclusione di un residuo e ogni scelta di mutazione deve avere una motivazione biofisica, strutturale o statistica documentata.*  
**Analisi Critica Completa dei Risultati:** Vedi **[`03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md)** e **[`SCALETTA_PROF/PUNTO_3_IPOTESI_MUTAZIONALI_E_SCREENING.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/SCALETTA_PROF/PUNTO_3_IPOTESI_MUTAZIONALI_E_SCREENING.md)**.

---

## 1. ✅ Stato Completato dello Screening Multi-Tool (`100%` Copertura dello Scaffold)

Tutti i 5 rami indipendenti del nostro **Orthogonal Consensus Funnel** sono stati completati e fusi in [`03_mutations/consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt):
1. **Dinamica Molecolare Vera del Wild-Type (`300 K` 1.5 ns vs `400 K` 2.5 ns):** identificati i **28 residui localmente instabili** ($\langle Q_i \rangle_{400\text{K}} < 0.60, N_0 \ge 30$, [`02_wt_baseline/unstable_residues_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/unstable_residues_report.txt)) e i legami a idrogeno che cedono ad alta temperatura ([`hbonds_details_300K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_300K.txt) vs [`hbonds_details_400K.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/hbonds_details_400K.txt)).
2. **FoldX 5 `PositionScan` sul `100%` dello Scaffold Modificabile (`88 su 88` posizioni $\times 20\text{ aa} = 1.760$ mutazioni):** completato grazie allo script parallelo a 8 core [`scripts/04b_run_remaining_foldx_parallel.sh`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/04b_run_remaining_foldx_parallel.sh) con scomposizione energetica a 15 termini (`03_mutations/energies_*_2MYO_Repair.txt`).
3. **Allineamento al Consenso Universale degli Ankyrin Repeats (*Mosavi et al., PNAS 2002*):** [`03_mutations/ankyrin_consensus_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ankyrin_consensus_report.txt) e [`06_presentation/fig6_ankyrin_consensus_alignment.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig6_ankyrin_consensus_alignment.pdf).
4. **FireProt 2.0 (Evolutivo MSA *Back-to-Consensus* `isBTC` + Energetico):** [`03_mutations/fireprot_results.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/fireprot_results.csv) (`2.242` righe).
5. **ThermoMPNN-D (Graph Neural Network su Singole e Doppie Epistatiche):** [`03_mutations/thermompnn_2MYO.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/thermompnn_2MYO.csv) (`2.242` singole) e [`03_mutations/thermompnn_double_2MYO.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/thermompnn_double_2MYO.csv) (`46.440` doppie).

---

## 2. ⚖️ Le 6 Regole Motivate di Selezione (Come Filtriamo da 118 Posizioni ai Candidati Finali)

Ogni scelta nel nostro imbuto decisionale segue **6 regole tassative e motivate**:

### Regola 1 — Esclusione dei 24 Residui della Functional Mask (`118` $\rightarrow$ `94` residui)
- **Residui esclusi:** `4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36, 41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100`.
- **Motivazione oggettiva:** Sono i 24 amminoacidi che nel complesso cristallografico [`3AAA`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/3AAA_complex_ref.pdb) hanno atomi pesanti a distanza **$< 4.5\text{ Å}$** dal bersaglio biologico (*Actin-Capping Protein*, inclusa la triade `TRP 8, ASP 44, GLU 78`, *Takeda et al., 2010*; [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) Slide 7 & 19). *(Questa regola ci ha permesso di scartare `D16W` e `I69G/M/R/N/S/H/K/L/D` proposte da FireProt 2.0).*

### Regola 2 — Esclusione delle Estremità N- e C-Terminali Libere (`94` $\rightarrow$ `88` residui interni)
- **Residui esclusi:** `1, 2, 3` (N-terminale) e `116, 117, 118` (C-terminale).
- **Motivazione oggettiva:** Le prime e ultime 3 posizioni sporgono liberamente nel solvente e fluttuano per natura termica senza contribuire alla barriera cooperativa di unfolding del solenoide Ankyrin (*Bekker et al., 2019* p. 433). *(Scarta `C2G`, `C2K`, `Q118W` di FireProt 2.0).*

### Regola 3 — Filtro di Sanità Biofisica (Scarto di Artefatti Chimici e Protezione dei Ponti H Inter-Repeat)
Anche se un singolo tool assegna $\Delta\Delta G < 0$, **scartiamo sempre**:
1. **Mutazioni che rompono la rete di legami a idrogeno inter-repeat delle Istidine (es. `S74A`):**  
   *Motivazione:* In Fase 1 ([`STATUS/discovery.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/STATUS/discovery.md#L15-L17)) abbiamo dimostrato che `SER 74 : OG` accetta il ponte a idrogeno inter-repeat ($2.44\text{ Å}$) da `HIS 40 (HSD)`. Mutare `SER 74 -> ALA` eliminerebbe l'ossigeno accettore `-OH` lasciando `HIS 40` desolvatata tra `ANK2` e `ANK3`!
2. **Mutazioni verso Cisteina (`-> CYS`, es. `T103C`, `T108C`):**  
   *Motivazione:* `2MYO` è priva di ponti disolfuro; introdurre un tiolo libero ($-\text{SH}$) espone a ossidazione e aggregazione.
3. **Mutazioni verso Prolina o Glicina (`-> PRO`, `-> GLY`) all'interno delle $\alpha$-eliche (es. `E52P` in `HELIX 4` o `T103G` nell'elica C-cap):**  
   *Motivazione:* Rompono il legame a idrogeno elicoidale $i \to i+4$ (*helix-breakers*, confermato da ThermoMPNN-D che boccia `E52P` con `+0.44` e `T103G` con `+0.69 kcal/mol`). `GLY` e `PRO` sono ammesse **solo nei loop/turn** (come **`H68G`**, **`A59G`**, **`D93G`**).

### Regola 4 — Motivazione Dinamica dalla MD del Wild-Type a `400 K` ($\langle Q_i \rangle_{400\text{K}}$ e Occupancy H-Bond)
- **Criterio:** La posizione scelta deve trovarsi in un residuo dinamicamente vulnerabile a `400 K` (**$\langle Q_i \rangle_{400\text{K}} < 0.60$** o nella fascia debole $< 0.72$ con alto RMSF / perdita di H-bond a `400 K`).

### Regola 5 — Convergenza Multi-Tool (5 Pilastri: MD 400K + FoldX + ThermoMPNN-D + FireProt 2.0 + Consenso Ankyrin)
- **Criterio:** Priorità assoluta alle mutazioni con **5/5 o 4/5 pilastri concordi** e meccanismo chimico-fisico chiaro nella scomposizione energetica di FoldX ([`fig8_consensus_energy_breakdown.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8_consensus_energy_breakdown.pdf)).

### Regola 6 — Complementarità Spaziale tra Repeat Distinti (`ANK2`, `ANK3`, `ANK4`)
- **Criterio:** Scegliere i mutanti singoli in moduli distinti del solenoide per massimizzare la stabilizzazione cooperativa a lungo raggio nel Doppio Mutante (*Lowe & Itzhaki, PNAS 2007*) ed evitare scontri sterici diretti.

---

## 3. 🏆 I Candidati Finalisti Selezionati e Validati in `FoldX BuildModel` (`n=5` Run, Slide 19)

| Variante | Modulo | Pilastri | $\Delta\Delta G_{\text{PositionScan}}$ | $\Delta\Delta G_{\text{BuildModel}}$ ($\mu \pm \sigma$) | $\langle Q_i \rangle_{400\text{K}}$ | ThermoMPNN-D | FireProt 2.0 | Sintesi del Meccanismo Fisico-Chimico & Epistasi |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **`K37L`** (`Mut1`) | **`ANK2`** | **5 / 5** | `-2.876` | **`-2.438 ± 0.870`** (`Best -3.32`) | **`0.4845` (UNSTABLE)** | **`-0.359`** | **`-0.89`** | Elimina la repulsione `+`/`+` con `ARG 36` e il ponte salino parassita `LYS 37 -> ASP 60` (`47.6%` a 300 K), stabilizzando il modulo più fragile a 400 K ($Q_{\text{ANK2}} = 0.586$). |
| **`T103L`** (`Mut2`) | **`ANK4`** | **5 / 5** | `-1.670` | **`-1.150 ± 0.143`** (`Best -1.25`) | **`0.5075` (UNSTABLE)** | **`-0.116`** | **`-1.18`** | Ripristina la `LEU` del consenso `TPLHLAA` all'ingresso dell'elica C-cap ($\Delta E_{\text{VdW+SolvH}} = -4.027\text{ kcal/mol}$), blindando il nucleo di folding primario C-terminale. |
| **`K37L + T103L`** (`Double`) | **`ANK2 + ANK4`** | **5/5 + 5/5** | `-4.546` (somma) | **`-3.652 ± 0.656`** (`Best -4.77`) | `0.485 & 0.508` | `-0.475` | `-2.07` | **100% Additiva in `BuildModel` ($\Delta\Delta\Delta G_{\text{coupling}} = -0.064\text{ kcal/mol}$!)**: rinforza contemporaneamente i due poli del solenoide che cedono a 400 K (`ANK2` e `ANK4`). |
| **`H68G`** (`Mut3 / BTC`) | **`ANK3`** | **4 / 5** | `-1.736` | **`-1.730 ± 0.009`** (`Best -1.75`) | **`0.7139` (WEAK)** | `+0.029` | **`-2.33` (`BTC, MAJ, RAT, HighC`)** | **#1 *Back-to-Consensus* di FireProt 2.0**: ripristina la `GLY` canonica in pos. 1 del $\beta$-turn `ANK2–ANK3`, rilassando l'outlier di Ramachandran sinistrorso ($\phi = +38.3^\circ$). |
| **`K37L + H68G + T103L`** (`Triple`) | **`ANK2+3+4`** | **5/5+4/5+5/5** | `-6.282` (somma) | **`-5.032 ± 0.210`** (`Best -5.39`) | `0.485, 0.714, 0.508` | `-0.446` | `-4.40` | **Record Assoluto (`-5.03 kcal/mol`):** combina 1 bullone per ciascuno dei 3 repeat (`ANK2`, `ANK3`, `ANK4`). |
