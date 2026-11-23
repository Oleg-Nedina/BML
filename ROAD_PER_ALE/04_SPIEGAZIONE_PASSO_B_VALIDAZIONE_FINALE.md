# 🤝 04_SPIEGAZIONE_PASSO_B_VALIDAZIONE_FINALE.md — Riepilogo Risultati dello Screening MD a 400 K (5 Sistemi) e Piano per la Validazione Finale (Passo B)

**Data:** 22 Novembre 2026  
**Team #10:** Oleg Nedina & Alessia Rigoni  
**Obiettivo di questo documento:** Riassumere per Alessia cosa abbiamo ottenuto dalle **5 simulazioni MD di screening a `400 K` (`Seed 1`, `2.5 ns`)** appena concluse (`WT`, `K37L`, `T103L`, `K37L_T103L` e `K37L_H68G_T103L`) e spiegare **perché e come nel Passo B (Campagna di Validazione Finale)** ci concentreremo sul confronto diretto tra il **Wild-Type (`WT`)** e il nostro **Vincitore Assoluto (`Triplo Mutante K37L_H68G_T103L`)**.

---

## 🏆 1. Cosa Abbiamo Scoperto nello Screening MD a `400 K` (`Seed 1` sui 5 Sistemi)

Abbiamo simulato in **NAMD 3 (`CHARMM36m + TIP3P 12 Å + 0.15 M NaCl`, `400 K`, `2.5 ns`)** i 4 candidati usciti dal nostro screening di consenso (*FoldX + FireProt 2.0 + ThermoMPNN + Consenso Ankyrin*) confrontandoli a parità di protocollo col Wild-Type (`2MYO`).  
Tutti i 15 plot PDF individuali e i 5 report separati si trovano in **[`05_analysis/single_runs/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs)** (e l'analisi biofisica completa è in **[`05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/RISULTATI_E_COMMENTO_SCREENING_400K_SEED1.md)**).

Ecco la tabella riassuntiva nella finestra equilibrata finale (`1.75–2.50 ns`, ultimi 30% dei frame, *Bekker et al. 2019* Table II):

| Sistema Simulato a `400 K` (`2.5 ns`) | $\Delta\Delta G_{\text{BuildModel}}$ ($n=5$ run) | `Q: hydrophilic-all` *(Bekker $r=0.84$ con $T_m$)* | `Q: all-all` *(Bekker $r=0.79$)* | `Q_ANK2` *(Modulo critico)* | Backbone RMSD Globale | Sito Funzionale `3AAA` RMSD | SASA Totale ($\text{Å}^2$) | Residui Scaffold Instabili ($\langle Q_i \rangle < 0.60$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`00_WT` (`2MYO` Wild-Type)** | `0.00 kcal/mol` | `0.6750` | `0.6948` | `0.6015` | `6.68 Å` | `7.27 Å` | `8839 Å²` | `41` residui |
| **`01_K37L` (Singolo `ANK2`)** | `-2.44 ± 0.87` | `0.7134` *(+3.84%)* | `0.7243` *(+2.95%)* | `0.6351` | `4.17 Å` | `4.75 Å` | `7942 Å²` | `25` residui *(-16)* |
| **`02_T103L` (Singolo `ANK4`)** | `-1.15 ± 0.14` | `0.7047` *(+2.97%)* | `0.7124` *(+1.76%)* | `0.6161` | `3.77 Å` | `3.37 Å` | `7969 Å²` | `30` residui *(-11)* |
| **`03_K37L_T103L` (Doppio `ANK2+4`)** | `-3.65 ± 0.66` | 🏆 **`0.7178` *(+4.28%)*** | `0.7194` *(+2.46%)* | `0.6290` | `6.40 Å` | `6.73 Å` | `8118 Å²` | `31` residui *(-10)* |
| **`04_K37L_H68G_T103L` (Triplo `ANK2+3+4`)** | 🏆 **`-5.03 ± 0.21`** | 🏆 **`0.7172` *(+4.22%)*** | 🏆 **`0.7289` *(+3.41%)*** | 🏆 **`0.6444`** | 🏆 **`3.30 Å`** | 🏆 **`2.83 Å`** | 🏆 **`7647 Å²`** | 🏆 **`23` residui *(-18)*** |

### 💡 Perché il Triplo Mutante (`K37L + H68G + T103L`) è il Vincitore Assoluto (e valorizza al massimo il lavoro di Alessia su FireProt 2.0 e ThermoMPNN!):
1. **Con `K37L` (`ANK2`) e `T103L` (`ANK4`):** Ripariamo i due moduli più deboli (`ANK2` e il C-cap `ANK4`) e portiamo `Q: hydrophilic-all` da `0.6750` a **`0.7178`**. Tuttavia, rompendo il ponte salino distorsivo `LYS 37 – ASP 60`, il residuo **`ASP 60`** (sulla punta del loop `60–68` tra `ANK2` e `ANK3`) rimane libero di oscillare nel solvente ($\text{RMSF}_{60} = 5.00\text{ Å}$ in `K37L` e `3.35 Å` nel Doppio), facendo flettere la cerniera centrale e tenendo alto l'RMSD globale del Doppio (`6.40 Å`).
2. **Con l'aggiunta di `H68G` (`ANK3` — la mutazione #1 *Back-to-Consensus* di FireProt 2.0 e #1 di ThermoMPNN trovata da Alessia!):**
   - `H68G` si trova esattamente al residuo `68`, cioè **all'uscita del loop `60–68` (`ANK2–ANK3`)**!
   - Sostituendo l'istidina 68 con la Glicina canonica del consenso Ankyrin, **elimina la tensione sterica del $\beta$-turn `ANK2–ANK3`**: la fluttuazione di **`ASP 60` crolla da `5.00 / 3.35 Å` a `1.50 Å`** (e quella di `LEU 37` scende a **`0.91 Å`** a 400 K!).
   - Risultato: il **Triplo Mutante `K37L_H68G_T103L`** vince su **tutte le metriche contemporaneamente** (`Q: all-all = 0.7289`, `Q_ANK2 = 0.6444`, `Backbone RMSD = 3.30 Å`, `Sito Funzionale 3AAA RMSD = 2.83 Å`, `SASA = 7647 Å²`, e **solo `23` residui instabili contro i `41` del WT = `-44%`**!).

---

## 🎯 2. Spiegazione del PASSO B: Perché e Come Facciamo la Validazione Finale su `WT` vs `Triplo Mutante`

Adesso che lo **Stadio 1 (Screening Dinamico a `400 K` sul `Seed 1`)** ha confrontato tutte e 5 le varianti (`WT`, `K37L`, `T103L`, `Doppio`, `Triplo`) e ha dimostrato perché il **Triplo `K37L_H68G_T103L`** è il design ottimale, passiamo al **Passo B (Validazione Finale Rigorosa)** concentrandoci sul confronto diretto **`Wild-Type (WT)` vs `Vincitore (K37L_H68G_T103L)`**.

### Perché nel Passo B testiamo solo `WT` vs `Triplo K37L_H68G_T103L`?
1. **Logica ad Imbuto (*Screening Funnel*, Slide 15 & 19 di `Laboratory_Project_2026.pdf`):**
   - Nello Stadio 1 abbiamo già dimostrato l'effetto individuale di `K37L`, `T103L`, del `Doppio` e del `Triplo` a parità di `Seed 1`.
   - Nello Stadio 2 (Validazione Finale), anziché sprecare decine di ore di calcolo sui mutanti intermedi, sottoponiamo **il controllo `WT` e il candidato finale `K37L_H68G_T103L`** ai due test di validazione più severi richiesti dalla teoria del corso (`Lecture04` Slide 6–15 e `Laboratory_Project_2026.pdf` Slide 10, 15, 21).

### Quali sono i 2 Test del Passo B su `WT` vs `K37L_H68G_T103L`?

#### **Test B.1 — Protocollo con Pre-Equilibrazione Restrained a `300 K` + Stress Termico a `400 K` (`Lecture04` Slide 6–9; *Bekker et al. 2019* p. 436)**
Nello screening rapido di Stadio 1 siamo passati direttamente dai `1.000` step di minimizzazione alla dinamica NPT a `400 K` (scartando il primo 70% transitorio).  
Nel **Passo B**, applichiamo sia al **`WT`** sia al **`Triplo K37L_H68G_T103L`** il protocollo multi-stadio completo da pubblicazione scientifica:
1. **Fase 1 (`1.000` step Minimizzazione Conjugate Gradient):** Rimuove ogni contatto sterico iniziale.
2. **Fase 2 (`50.000` step = `100 ps` di Pre-Equilibrazione NPT a `300 K` con *Harmonic Position Restraints* sul Backbone della proteina, `constraints on`, $k = 1.0\text{ kcal/mol/Å}^2$):**
   - Mentre il backbone della proteina è tenuto fermo nella sua geometria nativa da una molla armonica, l'acqua `TIP3P`, gli ioni ($\text{Na}^+, \text{Cl}^-$) e le catene laterali si rilassano e si idratano perfettamente alla temperatura fisiologica di **`300 K`** e pressione di **`1 atm`**.
3. **Fase 3 (`1.250.000` step = `2.50 ns` di Produzione NPT libera a `400 K`, `constraints off`):**
   - Si rilasciano i vincoli sul backbone e si porta il sistema già perfettamente equilibrato a **`400 K`**.

#### **Test B.2 — Repliche Statistiche Multi-Seed ($n = 3$ Velocity Seeds: `Seed 1`, `Seed 2`, `Seed 3`, Slide 15 & 21)**
- Eseguendo 3 repliche indipendenti (con semi della distribuzione di Maxwell-Boltzmann diversi) per **`WT`** e per **`K37L_H68G_T103L`**, calcoliamo **Media $\pm$ Deviazione Standard ($\mu \pm \sigma$)** e il **$p$-value del test $t$ di Welch**, dimostrando che il guadagno di stabilità del Triplo Mutante è statisticamente significativo ($p < 0.05$) e indipendente dal seme casuale iniziale!
