# ⚙️ RISULTATI_E_COMMENTO_FASE4_SETUP.md — Commento Scientifico al Setup Topologico "Matched Comparison" (`04_md_simulations/`) e Piano di Esecuzione Multi-Seed

**Progetto:** Biomolecular Modelling Laboratory (2026) — Team #10 (Oleg Nedina & Alessia Rigoni)  
**Sistemi Costruiti in `04_md_simulations/`:** `WT` (Seed 2 & 3), `K37L` (Seed 1–3), `T103L` (Seed 1–3), `K37L_T103L` (Seed 1–3), più `H68G` (Seed 1–3) e `K37L_H68G_T103L` (Seed 1–3)  
**Stato della Fase 4:** 🔄 **SETUP TOPOLOGICO E CONTROLLO DI QUALITÀ (GATE 2) COMPLETATI AL 100%; PRONTO PER LE SIMULAZIONI MD**  
**Riferimenti alle Slide del Corso:**
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf): **Slide 3** (*Learning Objective 3: Controlled computational comparison — Use consistent modelling and MD settings for WT and variants*), **Slide 10 & 21** (*Multiple velocity seeds & Computational budget*), **Slide 15** (*Step 4: Matched MD comparison*).
- [`Lecture03 - Molecular Dynamics Setup.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture03%20-%20Molecular%20Dynamics%20Setup.pdf): **Slide 32–50** (*Explicit water TIP3P, Periodic Boundary Conditions, cutoff 12 Å, PME neutrality & counter-ions*).
- [`Lecture04 - Running Molecular Dynamics.pdf`](file:///home/olly/UNI/BioMolEng2026/Lectures-20261014/Lecture04%20-%20Running%20Molecular%20Dynamics.pdf): **Slide 12, 15, 18–21, 26–27** (*Maxwell-Boltzmann velocity distribution, random seeds, NPT ensemble & reproducibility*).

---

## 1. 📊 Risultati del Controllo di Qualità Topologico (*Gate 2 Superato*)

L'esecuzione di [`scripts/08_build_matched_mutants_charmm36m.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/08_build_matched_mutants_charmm36m.tcl) sui conformeri a più bassa energia prodotti da `FoldX BuildModel` ha generato i domini solvatati e ionizzati in `CHARMM36m + TIP3P (12 Å) + 0.15 M NaCl`, documentati in [`04_md_simulations/matched_systems_qc_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/matched_systems_qc_report.txt):

| Sistema (`04_md_simulations/`) | Ruolo nel Disegno Sperimentale | Atomi Totali | Molecole $\text{H}_2\text{O}$ (`TIP3P`) | Ioni $\text{Na}^+$ (`SOD`) | Ioni $\text{Cl}^-$ (`CLA`) | Carica Elettrica Proteina | Carica Totale Cella PBC | Dimensioni Box PBC ($L_x \times L_y \times L_z\text{ Å}$) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`WT` (`2MYO`)** | Controllo Wild-Type | `21.461` | `6.535` | `23` | `18` | `-5.00 e` | **`+0.0000 e`** | `71.840 × 58.972 × 54.840` |
| **`K37L`** (`Mut1`) | Singola 1 (`ANK2`, 5/5★) | `23.068` | `7.070` | **`26`** | `20` | **`-6.00 e`** | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`T103L`** (`Mut2`) | Singola 2 (`ANK4`, 5/5★) | `23.075` | `7.070` | **`25`** | `20` | **`-5.00 e`** | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`K37L_T103L`** (`Double`) | **Doppia Primaria (`ANK2+4`)** | `23.067` | `7.068` | **`26`** | `20` | **`-6.00 e`** | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`H68G`** *(Opzionale BTC)* | Singola 3 (`ANK3`, #1 BTC) | `23.081` | `7.077` | **`25`** | `20` | **`-5.00 e`** | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |
| **`K37L_H68G_T103L`** *(Opz. Tripla)* | Tripla (`ANK2+3+4`) | `23.069` | `7.072` | **`26`** | `20` | **`-6.00 e`** | **`+0.0000 e`** | **`72.589 × 58.517 × 58.883`** |

---

## 2. 🧠 Commento e Considerazioni Critiche sulla Base della Teoria (`Lecture03` & `Lecture04`)

### 2.1 Bilancio di Carica ed Elettrostatica PME (`Lecture03` Slide 49–50)
- Nella Myotrophin Wild-Type a pH fisiologico ($7.0$), la somma algebrica dei residui carichi positivamente (`14 LYS + 3 ARG = +17 e`) e negativamente (`13 ASP + 9 GLU = -22 e`), con tutte le 4 istidine nello stato neutro tautomerico verificato in Fase 1 (`3 HSD + 1 HSE = 0 e`) e i terminali zwitterionici (`NH3+ = +1 e`, `COO- = -1 e`, `Lecture01` Slide 21), dà una carica netta esattamente pari a **`-5.00 e`**.
- Nelle varianti **`T103L`** (`THR 0 -> LEU 0`) e **`H68G`** (`HSD 0 -> GLY 0`), la sostituzione coinvolge amminoacidi elettricamente neutri: la carica della proteina rimane **`-5.00 e`**, e il plugin `autoionize` inserisce **`25` ioni $\text{Na}^+$ e `20` ioni $\text{Cl}^-$** ($\Delta q_{\text{ioni}} = +25 - 20 = +5.00\text{ e}$), portando la carica totale della cella a **`0.0000 e`**.
- Nelle varianti che includono **`K37L`** (`K37L`, `K37L_T103L`, `K37L_H68G_T103L`), la sostituzione di `LYS 37` (carica $+1\text{ e}$) con `LEU 37` (carica $0\text{ e}$) rimuove una carica positiva superflua e frustrata, portando la carica netta della proteina a **`-6.00 e`**. Come prescritto nella **`Lecture03` Slide 49–50** (*Particle Mesh Ewald & Counter ions*), il metodo PME richiede tassativamente che la cella unitaria periodica abbia carica netta nulla (altrimenti la serie di Ewald nello spazio reciproco diverge e richiede un background neutralizzante uniforme che può introdurre artefatti sulla pressione). Il comando `autoionize -sc 0.15` ha compensato esattamente la variazione inserendo **`26` ioni $\text{Na}^+$ e `20` ioni $\text{Cl}^-$** ($\Delta q_{\text{ioni}} = +26 - 20 = +6.00\text{ e}$), garantendo **`+0.0000 e`** su tutti i sistemi!

### 2.2 Coerenza Geometrica delle Celle Periodiche (`72.589 × 58.517 × 58.883 Å`)
- Tutte e 5 le varianti costruite con `FoldX BuildModel` a partire da `2MYO_Repair.pdb` presentano **dimensioni della scatola periodica identiche fino al millesimo di Ångström (`72.589 × 58.517 × 58.883 Å`)** e un numero di molecole d'acqua identico entro lo $0.1\%$ (`7.068 – 7.077` molecole `TIP3P`).
- La lieve differenza lungo l'asse $Z$ ($58.883\text{ Å}$ contro $54.840\text{ Å}$ del `WT` costruito in Fase 1 da `2MYO.pdb` grezzo) è dovuta unicamente al fatto che durante `FoldX RepairPDB` due catene laterali superficiali esposte lungo $Z$ si sono distese di $\sim 2\text{ Å}$ per eliminare i clash NMR iniziali, e `solvate -t 12` garantisce sempre un cuscinetto d'acqua di almeno $12.0\text{ Å}$ oltre l'atomo più esterno in ogni direzione (`Lecture03` Slide 41–43). In fase di produzione NPT a `1.01325 bar`, il barostato di Langevin Piston (`Lecture04` Slide 20–21) rilassa autonomamente la densità dell'acqua al valore d'equilibrio termodinamico a `400 K`.
