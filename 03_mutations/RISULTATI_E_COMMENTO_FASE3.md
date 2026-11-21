# 🎯 RISULTATI_E_COMMENTO_FASE3.md — Commento Scientifico, Screening Multi-Tool (`88/88` FoldX + FireProt 2.0 + ThermoMPNN-D), Risultati `BuildModel` e Considerazioni Teoriche della Fase 3 (`03_mutations/`)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Target:** Myotrophin (`2MYO_Repair.pdb`, 118 aa, `88` posizioni interne dello scaffold modificabili)  
**Stato della Fase 3:** ✅ **COMPLETATA E VALIDATA AL 100%** (inclusi `FoldX BuildModel` `n=5` con *matched WT references* ed estrazione dei PDB finali in `03_mutations/final_models_pdb/`)  
**Riferimenti alle Slide del Corso e Letteratura:**
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf): **Slide 3** (*Learning Objective 2: Mutation design — Explain how substitutions could improve local interactions*), **Slide 7** (*Computational filter on sequence space*), **Slide 15** (*Step 3: Mutation hypotheses and screening*), **Slide 16–19** (*FoldX energetic screening $\Delta\Delta G = \Delta G_{\text{mut}} - \Delta G_{\text{WT}}$; `RepairPDB` $\rightarrow$ `PositionScan` $\rightarrow$ `BuildModel` using matched WT references $\rightarrow$ `Selection for MD`*).
- [`Lecture05 - Beyond Conventional MD - Exploring and Designing Biomolecules.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture05%20-%20Beyond%20Conventional%20MD%20-%20Exploring%20and%20Designing%20Biomolecules.pdf): **Slide 12–16** (*Gibbs free energy & Entropy*), **Slide 48–59** (*Rosetta: rotamers, side-chain packing & Monte Carlo sampling*), **Slide 60–68** (*The protein design problem, Inverse Folding & ProteinMPNN*).
- *Mosavi LK, Minor DL, Peng ZY (2002).* **PNAS**, 99(25):16029–16034.
- *Lowe AR, Itzhaki LS (2007).* **PNAS**, 104(8):2679–2684.
- *Musil M et al. (2024).* **Nucleic Acids Research** (FireProt 2.0).
- *Dieckhaus H, Brocidiacono M, Randolph NZ, Kuhlman B (2024).* **PNAS**, 121(6):e2314853121 (ThermoMPNN / ThermoMPNN-D).

---

## 1. 📊 Sintesi dei Risultati Ottenuti nei 4 Step della Fase 3

### 1.1 Step 3.1 & 3.2 — Consenso Ankyrin e Riparazione Energetica `FoldX RepairPDB`
- **Allineamento ai 33 residui canonici degli Ankyrin Repeats ([`ankyrin_consensus_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ankyrin_consensus_report.txt) & [`fig6_ankyrin_consensus_alignment.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig6_ankyrin_consensus_alignment.pdf)):**
  - Ha rivelato che la Myotrophin naturale presenta specifiche deviazioni ("frustrazioni") rispetto al consenso termostabile `G-X-TPLHLAA-X-X-G-H-L-E-I-V-E-V-L-L-K-X-G-A-D-V-N-A-X-D` (*Mosavi et al., 2002*), spesso associate ad angoli di Ramachandran $(\phi, \psi)$ fuori regione nel `REMARK 500` di `2MYO.pdb`:
    - In **`ANK2`:** `LYS 37` al posto di `THR` in pos. 3 del motivo `TPLH`; `ALA 59`, `ASP 60` e `PRO 64` ($\psi = +121.5^\circ$) in tensione nel loop `GADVNA`.
    - In **`ANK3`:** `HIS 68` ($\phi = +38.3^\circ$) al posto della `GLY` universale in pos. 1 del $\beta$-turn; `LEU 73` al posto di `HIS` in pos. 6 (`TPLL`); `ASP 93` ($\phi = +35.4^\circ$), `LYS 94` (al posto di `VAL/ILE` in pos. 27) e `THR 95` (al posto di `ASN` in pos. 28) nel loop `GADKTV`.
    - In **`ANK4` (C-cap):** `THR 103` ($\phi = +8.6^\circ$) e `THR 108` ($\psi = -145.6^\circ$) all'ingresso e nel corpo dell'elica C-cap, oltre a `LEU 105` (`PHE 105` nell'uomo `3AAA`).
- **Effetto di `FoldX RepairPDB` ([`2MYO_Repair.fxout`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/2MYO_Repair.fxout)):**
  - L'ottimizzazione rotamerica di `RepairPDB` ha ridotto l'energia totale di FoldX da **`+319.97 kcal/mol`** (`2MYO.pdb` grezza) a **`+170.53 kcal/mol`** ([`2MYO_Repair.pdb`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/2MYO_Repair.pdb), **$\Delta E = -149.44\text{ kcal/mol}$**), abbattendo i clash di Van der Waals dell'80.4% (da `69.13` a `13.53 kcal/mol`) e i clash torsionali del 60.1% (da `86.95` a `34.67 kcal/mol`).

### 1.2 Step 3.3 & 3.4 — `FoldX PositionScan` (`88/88` Posizioni = `1.760` Varianti) + Unione con `FireProt 2.0` e `ThermoMPNN-D`
Incrociando i 5 pilastri indipendenti in [`03_mutations/consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt) (e analizzando tutte le 44 proposte di FireProt 2.0 in [`03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md)), sono emerse **due sole mutazioni singole in tutta la proteina con 5 Pilastri su 5 concordi e zero conflitti con la rete H-bond nativa**, più la **#1 mutazione evolutiva *Back-to-Consensus* di FireProt 2.0**:
1. **`K37L` in `ANK2` (5/5 Pilastri):** $\langle Q_{37} \rangle_{400\text{K}} = 0.4845$ (`UNSTABLE`), $\Delta\Delta G_{\text{PositionScan}} = -2.876\text{ kcal/mol}$, $\Delta\Delta G_{\text{ThermoMPNN}} = -0.359\text{ kcal/mol}$ (sito #1 assoluto dell'AI), $\Delta\Delta G_{\text{FireProt}} = -0.89\text{ kcal/mol}$.
2. **`T103L` in `ANK4` C-cap (5/5 Pilastri):** $\langle Q_{103} \rangle_{400\text{K}} = 0.5075$ (`UNSTABLE`), $\Delta\Delta G_{\text{PositionScan}} = -1.670\text{ kcal/mol}$ ($\Delta E_{\text{VdW+SolvH}} = -3.476\text{ kcal/mol}$), $\Delta\Delta G_{\text{ThermoMPNN}} = -0.116\text{ kcal/mol}$, $\Delta\Delta G_{\text{FireProt}} = -1.18\text{ kcal/mol}$, ripristina la `LEU` canonica del consenso `TPLHLAA`.
3. **`H68G` in `ANK3` $\beta$-turn (4/5 Pilastri, Rank #1 Score `-6.98` e #1 `BTC` FireProt):** $\langle Q_{68} \rangle_{400\text{K}} = 0.7139$ (`WEAK`, `RMSF = 3.63 Å`), $\Delta\Delta G_{\text{PositionScan}} = -1.736\text{ kcal/mol}$, $\Delta\Delta G_{\text{FireProt}} = -2.33\text{ kcal/mol}$ (`isBTC = true, majority = true, ratio = true, isHighCombined = true`), ripristina la `GLY` canonica in posizione 1 del repeat `ANK3` eliminando l'outlier di Ramachandran sinistrorso ($\phi = +38.3^\circ$).

### 1.3 Step 3.5 — Validazione Pre-MD con `FoldX BuildModel` (`n=5` Run con *Matched WT References*, Slide 19)
Eseguendo [`scripts/07_run_foldx_buildmodel.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/07_run_foldx_buildmodel.py) (`numberOfRuns=5` su 13 varianti singole, doppie e tripla contro i rispettivi *matched WT references* `WT_2MYO_Repair_*.pdb`, risultati in [`03_mutations/buildmodel_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/buildmodel_epistasis_report.txt) e [`06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf)), abbiamo ottenuto la prova regina per la scelta del Doppio Mutante:

| Variante (`BuildModel`, $n=5$) | Modulo | $\Delta\Delta G_{\text{BuildModel}}$ ($\mu \pm \sigma$, kcal/mol) | $\Delta\Delta G_{\text{Best}}$ (kcal/mol) | $\Delta E_{\text{VdW+SolvH}}$ | $\Delta E_{\text{Clash}}$ | $\Delta E_{\text{Entropy}}$ | Somma Singole ($\sum \Delta\Delta G_i$) | Accoppiamento ($\Delta\Delta\Delta G_{\text{coupling}}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`K37L`** (`Mut1`) | `ANK2` | **`-2.438 ± 0.870`** | **`-3.324`** | `-0.440` | `+0.030` | `-1.065` | — | — |
| **`T103L`** (`Mut2`) | `ANK4` | **`-1.150 ± 0.143`** | **`-1.250`** | **`-4.027`** | `+1.118` | `+1.375` | — | — |
| **`H68G`** (`Mut3 / BTC`) | `ANK3` | **`-1.730 ± 0.009`** | **`-1.747`** | `+0.555` | **`-0.598`** | **`-1.348`** | — | — |
| **`D93G`** | `ANK3` | `-1.616 ± 0.055` | `-1.674` | `+1.097` | `-1.048` | `-1.926` | — | — |
| **`T103F`** | `ANK4` | `-1.976 ± 0.319` | `-2.355` | `-4.589` | `+1.202` | `+1.467` | — | — |
| **`A59G`** | `ANK2` | `-2.555 ± 0.000` | `-2.555` | `+0.526` | `-1.722` | `-1.200` | — | — |
| **`D60L`** | `ANK2` | `-0.430 ± 1.743` | `-3.249` | `-3.446` | `+1.597` | `+0.182` | — | — |
| **`K37L + T103L`** (`Double A`) | **`ANK2 + ANK4`** | **`-3.652 ± 0.656`** | **`-4.766`** | **`-3.804`** | `+0.418` | `-0.272` | `-3.588` | **`-0.064` (100% Additiva!)** |
| **`H68G + T103L`** (`Double C`) | `ANK3 + ANK4` | **`-2.888 ± 0.418`** | **`-3.592`** | `-2.987` | `+0.688` | `-0.117` | `-2.880` | **`-0.008` (100% Additiva!)** |
| **`K37L + H68G`** (`Double B`) | `ANK2 + ANK3` | `-2.889 ± 0.076` | `-2.982` | `+0.001` | `-0.812` | `-1.810` | `-4.168` | `+1.278` (Sub-Additiva) |
| **`K37L + D93G`** (`Double D`) | `ANK2 + ANK3` | `-3.789 ± 0.858` | `-4.706` | `+1.487` | `-0.800` | `-3.530` | `-4.054` | `+0.264` (Quasi Additiva) |
| **`A59G + D60L`** (`Double E`) | `ANK2` locale | **`-3.548 ± 0.676`** | **`-4.236`** | **`-3.108`** | `+0.192` | `-0.478` | `-2.985` | **`-0.563` (Sinergia Epistatica!)** |
| **`K37L + H68G + T103L`** (`Triple`) | **`ANK2+3+4`** | **`-5.032 ± 0.210`** | **`-5.392`** | **`-2.968`** | **`-0.206`** | **`-1.596`** | `-5.318` | `+0.287` (Record `-5.03 kcal/mol`!) |

---

## 2. 🧠 Considerazioni Teoriche e Biofisiche: Cosa Ci Insegnano Questi Dati

### 2.1 Complementarità Termodinamica tra Entalpia ($\Delta H$) ed Entropia ($\Delta S$, `Lecture05` Slide 12–16)
Guardando la scomposizione energetica di `Average_2MYO_Repair.fxout` ([`buildmodel_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/buildmodel_epistasis_report.txt)), si vede come le nostre tre mutazioni vincitrici stabilizzino l'energia libera di Gibbs $\Delta G = \Delta H - T\Delta S$ attraverso **due meccanismi termodinamici ortogonali**:
1. **Stabilizzazione Entalpica di Core (`T103L` in `ANK4`):** Sostituire `THR 103` con `LEU 103` all'ingresso dell'elica C-cap paga un piccolo costo entropico di catena laterale ($\Delta E_{\text{Entropy}} = +1.375\text{ kcal/mol}$), ma guadagna un'enorme **energia entalpica di Van der Waals e solvatazione idrofobica ($\Delta E_{\text{VdW+SolvH}} = -4.027\text{ kcal/mol}$)**, sigillando il core C-terminale dall'ingresso dell'acqua!
2. **Stabilizzazione Entropica e di Rilascio di Strain (`H68G` in `ANK3` e `A59G` in `ANK2`):** Sostituire `HIS 68 -> GLY` (o `ALA 59 -> GLY`) in una posizione di $\beta$-turn con angolo diedro $\phi > 0$ rilassa i clash di backbone ($\Delta E_{\text{Clash}} = -0.598\text{ kcal/mol}$) e abbatte il costo entropico di mainchain (**$\Delta E_{\text{Entropy}} = -1.348\text{ kcal/mol}$**).
3. **Stabilizzazione Mista Elettrostatica/Solvatazione (`K37L` in `ANK2`):** Sostituire la Lisina carica frustrata `LYS 37` con `LEU 37` migliora simultaneamente la solvatazione ($\Delta E_{\text{SolvP+SolvH}} = -2.692\text{ kcal/mol}$) e l'entropia ($\Delta E_{\text{Entropy}} = -1.065\text{ kcal/mol}$).

### 2.2 Perché `K37L + T103L` è 100% Additiva (`-0.064 kcal/mol`) mentre `K37L + H68G` è Sub-Additiva (`+1.278 kcal/mol`)?
Questo risultato di `FoldX BuildModel` giustifica pienamente la **Regola 6 di Complementarità Spaziale** ([`MUTATION_SELECTION.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/MUTATION_SELECTION.md)):
- **`K37L` (in `ANK2`) e `H68G` (in `ANK3`)** distano solo $\sim 6\text{ Å}$ attraverso il ponte di **`HIS 67 (HSE)`** (che, come scoperto in Fase 1, forma un legame a idrogeno proprio con `GLY 35` e `LYS 37`!). Quando FoldX ottimizza `K37L` da sola o `H68G` da sola, entrambe rilassano lo stesso micro-intorno attorno a `HIS 67`; quando vengono combinate in `K37L + H68G`, quel rilassamento non può essere contato due volte ($\Delta\Delta\Delta G_{\text{coupling}} = +1.278\text{ kcal/mol}$).
- Al contrario, **`K37L` (nel polo N-terminale del core, `ANK2`) e `T103L` (nel polo C-terminale C-cap, `ANK4`)** agiscono su due moduli indipendenti e distanti ($d_{C_\alpha} > 20\text{ Å}$): i loro contributi termodinamici si sommano al 100% (**$\Delta\Delta G_{\text{K37L+T103L}} = -3.652 \pm 0.656\text{ kcal/mol}$**, `Best = -4.766 kcal/mol`, con accoppiamento $\Delta\Delta\Delta G_{\text{coupling}} = -0.064\text{ kcal/mol}$)!
- E quando nel **Triplo Mutante (`K37L + H68G + T103L`)** uniamo tutti e tre i repeat (`ANK2 + ANK3 + ANK4`), otteniamo una stabilizzazione complessiva di **`-5.032 ± 0.210 kcal/mol`** (`Best = -5.392 kcal/mol`)!

### 2.3 Conferma Fisica dell'Epistasi Predetta dall'Intelligenza Artificiale su `A59G + D60L`
Uno dei risultati scientificamente più eleganti dell'intero progetto è la convergenza tra **ThermoMPNN-D** (rete neurale Graph-based, `Lecture05` Slide 67) e **FoldX BuildModel** sulla coppia adiacente **`A59G + D60L`** ($d_{C_\alpha} = 3.82\text{ Å}$, all'ingresso del loop `GADVNA` di `ANK2`):
- In `BuildModel`, mutare **`D60L` da sola** dà un risultato incerto e debole ($\Delta\Delta G = -0.430 \pm 1.743\text{ kcal/mol}$ con clash $\Delta E_{\text{Clash}} = +1.597\text{ kcal/mol}$) perché il metile di `ALA 59` adiacente (outlier di Ramachandran $\phi = +175.1^\circ$) ostacola stericamente la catena laterale della Leucina 60.
- Quando però mutiamo contemporaneamente **`A59G + D60L`**, la sostituzione di consenso `A59G` libera lo spazio sterico nel backbone del turn, permettendo a `D60L` di impacchettarsi perfettamente ($\Delta E_{\text{VdW+SolvH}} = -3.108\text{ kcal/mol}$, clash ridotti a `+0.192 kcal/mol`) e sprigionando un **$\Delta\Delta G_{\text{BuildModel}} = -3.548 \pm 0.676\text{ kcal/mol}$** con una **sinergia epistatica positiva di $\Delta\Delta\Delta G_{\text{coupling}} = -0.563\text{ kcal/mol}$** — esattamente come predetto da ThermoMPNN-D (`epistasi = -1.292 kcal/mol` in [`consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt))!
