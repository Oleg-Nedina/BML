# 🧪 PASSO_B_VALIDAZIONE_FINALE.md — Protocollo Operativo per la Campagna di Validazione Finale (`WT` vs `Triplo K37L_H68G_T103L`)

**Data:** 22 Novembre 2026  
**Riferimenti Teorici:**
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 3** (*LO3: Controlled computational comparison; LO4: Assess uncertainty*), **Slide 10, 15, 21** (*Matched MD comparison with multiple velocity seeds*).
- [`Lecture04 - Running Molecular Dynamics.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture04%20-%20Running%20Molecular%20Dynamics.pdf) — **Slide 5–9** (*Energy minimization & Position-restrained pre-equilibration*), **Slide 12, 15, 18, 26–27** (*Maxwell-Boltzmann velocity seeds & statistical ensembles*).
- *Bekker, Ma & Kamiya (2019)*, **Protein Science**, 28:429–438 (Table II p. 434 & Methods p. 436).

---

## 🎯 1. Architettura a 2 Stadi della Validazione MD (Stadio 1 Completato $\rightarrow$ Stadio 2 / Passo B)

1. **Stadio 1 (Screening Dinamico a `400 K` sul `Seed 1` — ✅ COMPLETATO AL 100%):**
   - Confrontati tutti e 5 i sistemi (`00_WT`, `01_K37L`, `02_T103L`, `03_K37L_T103L`, `04_K37L_H68G_T103L`) a parità di protocollo e di `seed 12345` ([`05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md)).
   - Identificato il **Triplo Mutante `K37L_H68G_T103L`** (`ANK2 + ANK3 + ANK4`) come il **Vincitore Assoluto** su tutte le osservabili termodinamiche e strutturali (`ΔΔG = -5.03 kcal/mol`, `Q: all-all = 0.7289`, `Q: hydrophilic-all = 0.7172`, `Q_ANK2 = 0.6444`, `Backbone RMSD = 3.30 Å`, `Functional Mask 3AAA RMSD = 2.83 Å`, `23` residui instabili vs `41` del WT).
2. **Stadio 2 / Passo B (Validazione Finale Rigorosa Testa-a-Testa: `WT` vs `Triplo K37L_H68G_T103L`):**
   - Nel Passo B validiamo il confronto tra **`WT`** e **`K37L_H68G_T103L`** introducendo:
     1. **La Pre-Equilibrazione Restrained a `300 K` (`100 ps` = `50.000` step con `constraints on`, $k = 1.0\text{ kcal/mol/Å}^2$ sul backbone della proteina)** prima della fase di produzione non vincolata a **`400 K` (`2.50 ns`)**;
     2. **Le Repliche Multi-Seed ($n = 3$ velocity seeds indipendenti)** per calcolare **$\mu \pm \sigma$** e il **$p$-value statistico** tra `WT` e `K37L_H68G_T103L`.

---

## ⚙️ 2. Dettaglio Tecnico del Protocollo NAMD del Passo B (Pre-Equilibrazione Restrained a `300 K` + Produzione a `400 K`)

Per ciascuno dei due sistemi (`WT` e `K37L_H68G_T103L`), il protocollo con pre-equilibrazione esegue in un'unica chiamata NAMD 3 tre fasi fisiche consecutive:

1. **Creazione del file di riferimento dei vincoli armonici (`*_backbone_restraints.pdb`):**
   - Tutti gli atomi del solvente, degli ioni e delle catene laterali hanno colonna `beta = 0.00` (liberi di muoversi).
   - Gli atomi pesanti del backbone proteico (`name N CA C O`) hanno colonna `beta = 1.00` ($k = 1.0\text{ kcal/mol/Å}^2$).
2. **Stage 1 — Minimizzazione Energetica (`1.000` step con `constraints on`):**
   - Rilassa l'acqua `TIP3P`, gli ioni e le catene laterali attorno al backbone mantenuto nella conformazione nativa.
3. **Stage 2 — Pre-Equilibrazione NPT a Temperatura Fisiologica `300 K` (`50.000` step = `100 ps` con `constraints on`):**
   - `reinitvels 300`, `langevinTemp 300`, `langevinPistonTemp 300`, `run 50000`.
   - Il guscio di idratazione e la distribuzione ionica raggiungono l'equilibrio termico e barico a `300 K` senza alcuna distorsione precoce del backbone.
4. **Stage 3 — Produzione NPT a Stress Termico `400 K` (`1.250.000` step = `2.50 ns` con `constraints off`):**
   - Si azzera la costante di forza dei vincoli (`constraintScaling 0.0`), si porta il termostato e il barostato a **`400 K`** (`langevinTemp 400`, `langevinPistonTemp 400`, `reinitvels 400`) e si esegue la produzione di **`1.250.000` step (`2.50 ns`)**.
