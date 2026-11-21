# 🗺️ ROAD5.md — Fase 5: Analisi Statistica Multi-Seed (*Bekker Table II*), Ciclo Epistatico e Presentazione Finale (`05_analysis/` & `06_presentation/`)

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Target:** Confronto Statistico Multi-Seed a `400 K` (`WT` vs `K37L` vs `T103L` vs `K37L_T103L`)  
**Presentazione Finale:** 📅 **Mercoledì 16 Dicembre** (durata massima: **10 minuti**, Slide 4 di [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf))  
**Documenti Correlati:**
- [`ROADMAP/ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md) (Protocollo teorico dei 5 Livelli di Validazione e Go/No-Go Gates).
- [`SCALETTA_PROF/PUNTO_5_EVIDENZA_DECISIONE_E_PRESENTAZIONE.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/SCALETTA_PROF/PUNTO_5_EVIDENZA_DECISIONE_E_PRESENTAZIONE.md) (Monografia sul Punto 5 della scaletta del Prof e struttura delle 10 Slide).

---

## 📚 1. Fonti e Documentazione di Riferimento per la Fase 5

1. **Linee Guida Ufficiali del Progetto:**  
   [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 3** (*Learning Objective 4: Scientific interpretation — Assess uncertainty and discuss how computational screening can support experimental evaluation*), **Slide 4** (*Project outputs*), **Slide 11–13** (*Q-value correlation with $T_m$*), **Slide 15** (*Step 5: Evidence and decision. Compare Q and supporting observables. If necessary re-iterate*), **Slide 22–23** (*RMSF via `analysisMD.tcl` & Q-value via `VMDextensions.tcl`*).
2. **Paper di Riferimento:**  
   *Bekker GJ, Ma B, Kamiya N (2019).* **Protein Science**, 28(2):429–438 — **Table II (p. 434)** e **Tables S1, S8–S10 (pp. 433–434)**.
3. **Letteratura sul Folding Cooperativo della Myotrophin (`2MYO`):**  
   *Lowe AR, Itzhaki LS (2007).* **PNAS**, 104(8):2679–2684.

---

## 🛠️ 2. Pipeline Automatica di Analisi e Validazione (Pronta in `scripts/`)

Non appena termina anche solo il primo batch della Fase 4 (`seed1` o tutte le 3 repliche), l'intera analisi della Fase 5 si esegue con **2 soli comandi** dalla cartella `BML/`:

### 💻 Comando 5.1 — Estrazione delle Osservabili in VMD ([`scripts/09_analyze_matched_validation.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/09_analyze_matched_validation.tcl))
```bash
vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl
```
- **Cosa calcola automaticamente per ogni replica `.dcd` presente su disco:**
  1. Ricomposizione PBC (`pbc wrap`) e allineamento sul backbone del frame 0 (`measure fit`, Slide 22).
  2. **Backbone RMSD($t$)** e **RMSD ristretto ai 24 residui della Functional Mask (`3AAA`)** (per dimostrare che il sito di legame per l'Actin-Capping Protein resta intatto!).
  3. **RMSF dei $C_\alpha$** per ciascuno dei 118 residui.
  4. **Hydrogen Bonds** totali nel tempo + dettaglio di occupancy (`-type unique`, `05_analysis/*_hbonds_details.txt`).
  5. **SASA($t$)**, **Raggio di Girazione $R_g(t)$** e **% di $\alpha$-elica STRIDE** ([`analysisMD.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/analysisMD.tcl)).
  6. **$Q$-value nel tempo** per le **3 classi chimiche di Bekker et al.** (`Q_all`, `Q_phil_all` [$r=0.84$ con $T_m$!], `Q_phob_small`), per i **4 singoli Ankyrin Repeats** (`Q_ANK1..4`) e per le **3 interfacce inter-repeat** (`Q_int12`, `Q_int23`, `Q_int34`).
  7. **$\langle Q_i \rangle_{400\text{K}}$ per singolo residuo (`1..118`)** sugli ultimi 30% dei frame (`876..1251`) e conteggio dei residui dello scaffold localmente instabili ($\langle Q_i \rangle < 0.60, N_0 \ge 30$).

### 💻 Comando 5.2 — Compilazione della *Bekker Table II*, Test Statistici e Grafici Vettoriali `.pdf` ([`scripts/10_plot_final_validation_bekker.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/10_plot_final_validation_bekker.py))
```bash
python3 scripts/10_plot_final_validation_bekker.py
cat 05_analysis/bekker_table2_and_epistasis_report.txt
```
- **Output generati:**
  - **`05_analysis/bekker_table2_and_epistasis_report.txt`**:
    1. **Table II (stile *Bekker et al., 2019* p. 434):** Media $\pm$ SD ($\mu \pm \sigma$) tra i seed di `RMSD_BB`, `Q: hydrophilic-all`, `Q: all-all`, `Q: hydrophobic-small`, numero di residui instabili dello scaffold ($Q_i < 0.60$) e $p$-value del test $t$ di Welch rispetto al WT.
    2. **Tabella Repeat-Resolved & `analysisMD.tcl`:** Media $\pm$ SD di `Q_ANK1..4`, `Q_int(1-2, 2-3, 3-4)`, `% Alpha-Helix`, $R_g$ e `RMSD della Functional Mask 3AAA`.
    3. **Ciclo Epistatico Dinamico del Doppio Mutante:** confronto tra $\Delta Q_{\text{K37L+T103L}}$ e $\Delta Q_{\text{K37L}} + \Delta Q_{\text{T103L}}$.
  - **Grafici Vettoriali `.pdf` in `06_presentation/`:**
    - `06_presentation/fig9_matched_q_time_and_repeats.pdf`
    - `06_presentation/fig10_per_residue_delta_qi_and_rmsf.pdf`
    - `06_presentation/fig11_epistatic_cycle_and_observables.pdf`

---

## 📊 3. Riepilogo Completo delle Figure Vettoriali `.pdf` in `06_presentation/` per le Slide del 16 Dicembre

| Figura in `06_presentation/` | Stato Attuale | Punto della Scaletta Prof (`Slide 15`) | Cosa Mostra nelle Slide |
| :--- | :---: | :---: | :--- |
| **[`fig5_wt_baseline_plots.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5_wt_baseline_plots.pdf)** | ✅ **Pronta (Dati Reali)** | **Punto 2** (Slide 4) | 4 quadranti WT `300 K` vs `400 K`: Backbone RMSD, $Q(t)$, RMSF($C_\alpha$) e profilo $\langle Q_i \rangle_{400\text{K}}$ con i 28 residui instabili. |
| **[`fig5b_wt_analysisMD_observables.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig5b_wt_analysisMD_observables.pdf)** | ✅ **Pronta (Dati Reali)** | **Punto 2** (Slide 4) | 4 quadranti `analysisMD.tcl` WT: `% Alpha-Helix` STRIDE, Hydrogen Bonds, Raggio di Girazione $R_g$ e $Q(t)$ per `ANK1–4`. |
| **[`fig6_ankyrin_consensus_alignment.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig6_ankyrin_consensus_alignment.pdf)** | ✅ **Pronta (Definitiva)** | **Punto 1 & 3** (Slide 5) | Matrice dei 4 Ankyrin Repeats allineati ai 33 residui del consenso di *Mosavi et al. (2002)* con Functional Mask `3AAA`. |
| **[`fig7_foldx_position_scan_heatmap.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig7_foldx_position_scan_heatmap.pdf)** | ✅ **Pronta (Copertura 88/88)** | **Punto 3** (Slide 6) | Classifica delle Top 20 mutazioni stabilizzanti da FoldX `PositionScan` sull'intero scaffold modificabile. |
| **[`fig8_consensus_energy_breakdown.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8_consensus_energy_breakdown.pdf)** | ✅ **Pronta (5 Pilastri)** | **Punto 3** (Slide 6) | Scomposizione fisico-chimica ($\Delta E_{\text{VdW+SolvH}}$, $\Delta E_{\text{Clash}}$, $\Delta E_{\text{HBond}}$) dei migliori candidati multi-tool. |
| **[`fig8b_buildmodel_singles_doubles_epistasis.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf)** | ✅ **Pronta (`n=5` Run)** | **Punto 3** (Slide 6) | $\Delta\Delta G_{\text{BuildModel}}$ ($\mu \pm \sigma$ contro *matched WT references*, Slide 19) e ciclo di additività/epistasi statica. |
| **`fig9_matched_q_time_and_repeats.pdf`** | ⏳ *Generata da Script 10* | **Punto 4 & 5** (Slide 7) | Curve $Q(t)$ multi-seed ($\mu \pm \sigma$) e barplot comparativi per classi chimiche e repeat (`WT` vs `K37L`, `T103L`, `K37L+T103L`). |
| **`fig10_per_residue_delta_qi_and_rmsf.pdf`** | ⏳ *Generata da Script 10* | **Punto 5** (Slide 8) | Guarigione dei residui instabili ($\Delta\langle Q_i \rangle_{400\text{K}}$) e riduzione della flessibilità ($\Delta\text{RMSF}_{C_\alpha}$) lungo i 118 residui. |
| **`fig11_epistatic_cycle_and_observables.pdf`** | ⏳ *Generata da Script 10* | **Punto 5** (Slide 9) | Ciclo epistatico dinamico del doppio mutante e osservabili `analysisMD.tcl` (`% Alpha-Helix`, $R_g$, `Functional Mask RMSD`). |
