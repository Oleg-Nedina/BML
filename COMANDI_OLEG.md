# 🛠️ COMANDI_OLEG.md — Comandi Operativi Mutazioni (3 Seeds a 400 K), Generazione Risultati e Consegna ad Alessia

**Corso:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Ruolo Oleg in questa fase:** Eseguire il triplicato a `400 K` (`3 velocity seeds`: `12345, 23456, 34567`, `2.5 ns`) sulle mutazioni, estrarre tutte le metriche (`Q`, `RMSD`, `RMSF`, `SASA`, `Rgyr`, `H-Bonds 3.5 Å / 30°`, `7 α-eliche H1–H7` con varianza $\mu \pm \sigma$), analizzare i risultati rispetto alla baseline `WT` e inoltrare i pacchetti pronti ad Alessia.

---

## 🚀 1. Comandi per Lanciare il Triplicato (`3 Seeds` a `400 K`) di Ogni Mutazione

Esegui dalla cartella `/home/olly/UNI/BioMolEng2026/laboratorio/BML`.  
*(Tutti e 4 gli script usano `namd3 +p8`, saltano automaticamente il **`Seed 1` (`12345`) già completato** ed eseguono **`Seed 2` (`23456`)** e **`Seed 3` (`34567`)**).*

### 🏆 Priorità 1 — Triplo Mutante (`K37L_H68G_T103L`, vincitore assoluto)
```bash
bash 04_md_simulations/K37L_H68G_T103L/run_K37L_H68G_T103L_3seeds.sh 8
```

### 🥈 Priorità 2 — Doppio Mutante (`K37L_T103L`)
```bash
bash 04_md_simulations/K37L_T103L/run_K37L_T103L_3seeds.sh 8
```

### 🥉 Priorità 3 — Singola Mutazione 1 (`K37L` in Repeat 1)
```bash
bash 04_md_simulations/K37L/run_K37L_3seeds.sh 8
```

### 4️⃣ Priorità 4 — Singola Mutazione 2 (`T103L` in C-terminal)
```bash
bash 04_md_simulations/T103L/run_T103L_3seeds.sh 8
```

*(Per monitorare l'avanzamento in un secondo terminale, es. per il Triplo Seed 2:)*
```bash
grep -E "TIMING|Benchmark|End of program|FATAL ERROR" 04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed2.log | tail -n 10
```

---

## 📊 2. Comandi per Analisi Completa e Generazione dei Risultati da Mandare ad Alessia

Appena finisce una mutazione (o dopo aver fatto `git pull` per ricevere i nuovi seed del `WT` di Alessia), lancia questa sequenza dalla cartella `BML/`:

### A. Estrazione Dati VMD (`Q-value` *Yang et al. 1998* + `Frame 0 = DCD`, `RMSD`, `RMSF`, `SASA`, `Rgyr`, `H-Bonds 3.5 Å / 30°`, `7 α-eliche H1–H7`):
```bash
vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl && \
vmd -dispdev text -e scripts/12_analyze_hbonds_and_helices.tcl
```
*(Entrambi hanno cache intelligente: analizzano solo i `.dcd` nuovi o aggiornati!).*

### B. Generazione Tabelle ($\mu \pm \sigma$ e $p$-value) e Tutti i Grafici Comparativi (`.pdf` + `.png`):
```bash
python3 scripts/10_plot_single_runs_400K.py && \
python3 scripts/10_plot_final_validation_bekker.py && \
python3 scripts/11_compare_wt_vs_triple_for_prof.py && \
python3 scripts/12_plot_hbonds_and_helices.py
```

### C. Stampa a schermo delle 3 Tabelle Riepilogative da verificare e inoltrare ad Alessia:
```bash
cat 05_analysis/bekker_table2_and_epistasis_report.txt
cat 06_presentation/presentazione_prof_WT_vs_Triplo/00_TABELLA_CONFRONTO_WT_vs_TRIPLO.txt
cat 06_presentation/presentazione_prof_WT_vs_Triplo/05_hbonds_and_helices_summary.txt
```

### D. Visualizzazione 3D in VMD (per eventuali screenshot extra per le slide):
```bash
vmd -e scripts/viz_04_mutations_zoom.tcl          # Zoom 3D statico WT vs K37L, H68G, T103L (Yang 1998)
vmd -e scripts/viz_05_md_comparison_400K.tcl      # Animazione MD 400 K del Triplo (o -args <VAR> <seed>)
vmd -e scripts/viz_06_wt_vs_triple_400K.tcl       # Sovrapposizione 3D WT (Rosso) vs Triplo (Viola) + H6/H7
```

### E. Mappa Ordinata delle Cartelle di Output (dove salva ogni script)
- `04_md_simulations/<VAR>/` $\rightarrow$ topologia (`.psf`/`.pdb`), config, `.dcd`, `.log` e `run_<VAR>_3seeds.sh`
- `05_analysis/` $\rightarrow$ dati grezzi `.dat` per replica (`<VAR>_s<1..3>_*`) + `bekker_table2_and_epistasis_report.txt`
- `05_analysis/single_runs/<00_WT .. 04_K37L_H68G_T103L>/` $\rightarrow$ 3 PDF + 1 report dedicati per ogni singola variante
- `05_analysis/helices_hbonds/` $\rightarrow$ dati H-bonds (`3.5 Å / 30°`), matrici STRIDE e geometria `H1–H7`
- `06_presentation/` $\rightarrow$ figure master comparate (`fig5`..`fig11`)
- `06_presentation/presentazione_prof_WT_vs_Triplo/` $\rightarrow$ pacchetto pronto per le slide (`00_TABELLA...`, `01`–`03`, `05`–`06` in `.pdf` e `.png`)

---

## 📋 3. Scaletta Operativa di Oleg (Workflow Passo-Passo)

1. **Lancio MD Mutazioni (Triplicato 3 Seeds):** Partire dal Triplo `K37L_H68G_T103L`, poi Doppio `K37L_T103L`, poi Singoli `K37L` e `T103L`.
2. **Check Qualità dei Log:** Verificare `End of program` e assenza di `FATAL ERROR` su ogni `.log` generato.
3. **Estrazione & Confronto Multi-Metrica contro Baseline WT:**
   - Lanciare i comandi della Sezione 2 (`09` + `12` VMD $\rightarrow$ `10` + `11` + `12` Python).
   - Verificare **media $\pm$ deviazione standard ($\mu \pm \sigma$)** e **$p$-value** tra i 3 seed su:
     - $Q$-value (`all-all`, `hydrophilic-all`, `hydrophobic-small`, 4 regioni *Yang et al. 1998*);
     - `Backbone RMSD`, `Sito Funzionale 3AAA RMSD`, `Rgyr`, `SASA`, `RMSF`;
     - **Nuove metriche Prof.:** `H-Bonds (3.5 Å / 30°)` totali e backbone + **comportamento delle 7 $\alpha$-eliche (`H1–H7`, con focus su `H3–H4` e le ultime eliche `H6–H7`)**.
4. **Sintesi & Consegna ad Alessia (`git push`):**
   - Aggiornare le considerazioni nei report (possiamo farlo insieme appena escono i numeri del triplicato!) e pushare tutto su GitHub così Alessia trova già pronti i grafici e le tabelle in `06_presentation/` e `06_presentation/presentazione_prof_WT_vs_Triplo/`:
     ```bash
     git add 05_analysis/ 06_presentation/ *.md ROADMAP/ ROAD_PER_ALE/
     git commit -m "Updated 3-seed MD analysis, H-bonds and alpha-helices H1-H7 comparison"
     git pull --no-rebase --no-edit && git push
     ```
