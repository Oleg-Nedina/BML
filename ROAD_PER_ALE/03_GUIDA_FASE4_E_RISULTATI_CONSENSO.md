# 🏆 03_GUIDA_FASE4_E_RISULTATI_CONSENSO.md — Cosa Abbiamo Trovato Unendo i Tuoi CSV (FireProt 2.0 + ThermoMPNN-D) con la MD a 400 K e Come Dividiamo la Fase 4

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Data:** 21 Novembre 2026  

---

## 🎉 1. I Risultati Straordinari dell'Unione dei Tuoi File `.csv` con i Calcoli di Oleg

I 3 file che hai generato e caricato in `03_mutations/` ([`fireprot_results.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/fireprot_results.csv), [`thermompnn_2MYO.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/thermompnn_2MYO.csv) e [`thermompnn_double_2MYO.csv`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/thermompnn_double_2MYO.csv)) sono stati uniti con:
- La **Dinamica Molecolare del Wild-Type a `400 K`** (`2.5 ns`, [`02_wt_baseline/unstable_residues_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/02_wt_baseline/unstable_residues_report.txt)),
- Lo **screening FoldX su tutte le `88/88` posizioni modificabili** ([`03_mutations/consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt)),
- E la costruzione 3D **`FoldX BuildModel` (`n=5` run indipendenti)** ([`03_mutations/buildmodel_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/buildmodel_epistasis_report.txt)).

Ecco le **4 scoperte principali** che porterai nelle tue Slide il 16 Dicembre (descritte in modo esteso in [`03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md)):

1. **Abbiamo trovato 2 mutazioni singole con `5 Pilastri su 5` concordi in due moduli opposti della proteina (`ANK2` e `ANK4`):**
   - **🥇 `K37L` (nel Repeat 2, `ANK2`):**
     - Nella MD a `400 K`, `ANK2` è il modulo che cede di più ($Q_{\text{ANK2}} = 0.586$) e `LYS 37` perde oltre metà dei suoi contatti ($\langle Q_{37} \rangle = 0.4845$).
     - Nel tuo **ThermoMPNN-D**, `K37` è **il sito #1 più stabilizzante di tutta la proteina** (`K37L = -0.359 kcal/mol`)!
     - Nel tuo **FireProt 2.0** dà **`-0.89 kcal/mol`**, in **FoldX `PositionScan`** dà **`-2.876 kcal/mol`** e in **`FoldX BuildModel` (`n=5`)** dà **`-2.438 ± 0.870 kcal/mol`** (`Best = -3.324 kcal/mol`).
   - **🥈 `T103L` (nel C-cap, `ANK4`):**
     - Nella MD a `400 K`, l'ingresso dell'elica del C-cap (`LEU 102 – THR 103`, che secondo *Lowe & Itzhaki, PNAS 2007* è il **nucleo di folding primario** della Myotrophin) si srotola ($\langle Q_{103} \rangle = 0.5075$).
     - `T103L` ripristina la **Leucina (`L`) canonica del Consenso Ankyrin** (*Mosavi et al., 2002*), è approvata dal tuo **FireProt 2.0 (`-1.18 kcal/mol`)**, dal tuo **ThermoMPNN-D (`-0.116 kcal/mol`)**, da **FoldX `PositionScan` (`-1.670 kcal/mol`)** e in **`FoldX BuildModel` (`n=5`)** dà **`-1.150 ± 0.143 kcal/mol`** con un guadagno idrofobico nel core di ben **`-4.027 kcal/mol`**!
2. **Il Doppio Mutante `K37L + T103L` (`ANK2 + ANK4`) è `100% Additivo` in `FoldX BuildModel`:**
   - In `BuildModel` (`n=5`), la somma teorica delle due singole è `-2.438 + (-1.150) = -3.588 kcal/mol`, e il Doppio Mutante **`K37L + T103L`** ottiene **`-3.652 ± 0.656 kcal/mol`** (`Best = -4.766 kcal/mol`), con un accoppiamento epistatico $\Delta\Delta\Delta G_{\text{coupling}} = \mathbf{-0.064\text{ kcal/mol}}$ (additività perfetta e zero interferenza sterica!).
3. **La tua #1 *Back-to-Consensus* di FireProt 2.0 (`H68G` in `ANK3`) e il Triplo Mutante da Record (`K37L + H68G + T103L`):**
   - Nel tuo `fireprot_results.csv`, **`H68G`** è l'unica mutazione con tutti e 4 i flag evolutivi ed energetici attivi (`isBTC = true`, `majority = true`, `ratio = true`, `isHighCombined = true`, $\Delta\Delta G_{\text{FireProt}} = -2.33\text{ kcal/mol}$, $\Delta\Delta G_{\text{BuildModel}} = -1.730 \pm 0.009\text{ kcal/mol}$) perché ripristina la Glicina canonica in posizione 1 del $\beta$-turn `ANK2–ANK3`.
   - Unendo `K37L` (`ANK2`) + `H68G` (`ANK3`) + `T103L` (`ANK4`), il **Triplo Mutante** raggiunge in `BuildModel` il record assoluto di **`-5.032 ± 0.210 kcal/mol`** (`Best = -5.392 kcal/mol`)!
4. **Due "Chicche" Scientifiche d'Oro per le Tue Slide (Slide 5 & 6):**
   - **Perché abbiamo scartato `S74A` (che piaceva a FireProt, FoldX e ThermoMPNN)?** Perché in Fase 1 avevamo scoperto che l'ossigeno `OG` di `SER 74` accetta il ponte a idrogeno inter-repeat ($2.44\text{ Å}$) da `HIS 40 (HSD)`! Mutare `SER 74 -> ALA` toglierebbe il gruppo `-OH` rompendo la cucitura tra `ANK2` e `ANK3`.
   - **Conferma fisica dell'epistasi AI su `A59G + D60L`:** Il tuo file `thermompnn_double_2MYO.csv` aveva previsto una forte sinergia epistatica locale per la coppia adiacente `A59G + D60L` in `ANK2` (`-0.756 kcal/mol` doppia contro `+0.535 kcal/mol` somma delle singole). Quando l'abbiamo testata in `FoldX BuildModel`, anche FoldX ha confermato una **sinergia epistatica positiva di `-0.563 kcal/mol`** (`-3.548 kcal/mol` doppia contro `-2.985 kcal/mol` somma delle singole), perché `A59G` sblocca lo spazio sterico per `D60L`!

---

## 🖥️ 2. Come Funziona Adesso la Fase 4 (Simulazioni MD Multi-Seed a `400 K`)

Oleg ha già generato e pushato nella cartella **`04_md_simulations/`** tutte le scatole d'acqua `CHARMM36m + TIP3P + 0.15 M NaCl` e i file `.conf` pronti per NAMD 3 per i 3 velocity seeds (`seed1`, `seed2`, `seed3` a `400 K`, `2.5 ns` ciascuno).

Se volete dividere a metà le simulazioni dei `Seed 2` e `Seed 3` tra il PC Linux di Oleg e il tuo PC Windows:
1. Fai `git pull` nella cartella `BML`.
2. Sul tuo PC Windows (da Git Bash o PowerShell dentro `BML/`), puoi lanciare le 4 simulazioni assegnate a `alessia` (`WT seed 2 & 3` + `T103L seed 2 & 3`):
   ```bash
   bash 04_md_simulations/run_batch_phase4.sh alessia 8
   ```
   *(Sostituisci `8` con il numero di core fisici della tua CPU Windows; lo script salta da solo le run già fatte!).*
