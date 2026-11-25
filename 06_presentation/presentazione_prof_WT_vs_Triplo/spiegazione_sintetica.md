# 🎙️ SPIEGAZIONE SINTETICA PASSO-PASSO DEL PROGETTO (Traccia Completa per il Colloquio con il Prof. Gautieri)

**Corso:** Biomolecular Modelling Laboratory (2026) — Prof. Alfonso Gautieri  
**Team #10:** Oleg Nedina & Alessia Rigoni  
**Proteina Assegnata:** Myotrophin (*Rattus norvegicus*, PDB ID: [`2MYO`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/2MYO.pdb), 118 amminoacidi — 4 moduli Ankyrin-like)

---

## 📌 PARTE 1 : Analisi Strutturale, Funzione Biologica e *Functional Mask* (`2MYO` vs `3AAA`)

### 🗣️ Parte 1 : Discorso da dire al Professore
> *"Professore, siamo partiti dall'analisi strutturale e biologica della nostra proteina, la **Myotrophin (`2MYO`)** di ratto, una piccola proteina intracellulare di **118 amminoacidi** (~13 kDa) appartenente alla famiglia delle proteine a **Ankyrin Repeat**.*
> *A livello biologico, la sua funzione meglio caratterizzata sperimentalmente è la regolazione del citoscheletro di actina: la Myotrophin si lega direttamente all'**Actin-Capping Protein (CP)** e ne inibisce per ingombro sterico l'attività di capping sulle estremità dei filamenti di actina. Di conseguenza, il nostro obiettivo di Protein Engineering è **aumentare la stabilità termica dello scaffold Ankyrin senza alterare la superficie di riconoscimento molecolare per CP**.*
> *La struttura depositata in `2MYO` è un modello **Solution NMR** (struttura media minimizzata da 90 conformeri a 297 K e pH 6.2), monomerico (`Chain A`), composto da **7 segmenti $\alpha$-elicoidali** (~46.6% della sequenza secondo STRIDE) organizzati in **4 regioni modulari**: una regione N-terminale (`residui 1–27/34`), **due Ankyrin repeats centrali completi** (`repeat 1: 28–60 / 35–67` e `repeat 2: 61–93 / 68–100`, ciascuno con il classico motivo elica–turn–elica e un loop $\beta$-hairpin sporgente alle posizioni `32–36` e `65–69`), e una regione C-terminale (`94–118`).*
> *Per definire in modo oggettivo quali residui non vanno assolutamente mutati, abbiamo sovrapposto `2MYO` alla struttura cristallografica del complesso Myotrophin–Actin Capping Protein (**PDB `3AAA`**, `Chain C`, $\text{RMSD} = 2.365\text{ Å}$, identità di sequenza 99.2%) e abbiamo inserito in una **Functional Mask intoccabile tutti i 24 residui a distanza $< 4.5\text{ Å}$ dal partner di legame CP**."*

### ✅ Perché è valido e in che modo oggettivo l'abbiamo checkato
- **Check di completezza atomica ([`report/01_scripts/01_check_missing_residues.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/01_scripts/01_check_missing_residues.tcl) & [`scripts/01_inspect_2myo.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/01_inspect_2myo.tcl)):** Confrontati i record `SEQRES` e `ATOM`: `118/118` residui continui (`MET 1` $\rightarrow$ `GLN 118`), `118` atomi $C_\alpha$, `902` atomi pesanti (`0` atomi mancanti, `0` duplicati, `0` altLoc) e `913` idrogeni NMR (`1.815` atomi totali).
- **Check ponti disolfuro:** Misurate le distanze tra gli atomi di zolfo `SG` delle 3 cisteine (`CYS 2, CYS 45, CYS 83`): la distanza minima è **`9.434 Å`** (`CYS 45 – CYS 83`), ben oltre i `~2.05 Å` di un ponte disolfuro $\rightarrow$ **tutte e 3 le cisteine sono nello stato ridotto (`-SH`)**.
- **Check impacchettamento idrofobico statico (SASA/RSA, [`report/01_scripts/02_sasa_per_residue.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/01_scripts/02_sasa_per_residue.tcl)):** Calcolata la Relative Solvent Accessibility (normalizzata secondo *Tien et al. 2013*), trovando **19 residui completamente sepolti (`RSA < 0.05`)** concentrati nel core idrofobico dei repeat centrali (es. `ALA 42, ILE 50, LEU 51, LEU 54`), más **21 legami a idrogeno** intramolecolari e **2 ponti salini statici** entro $4.0\text{ Å}$ (`ASP 27 – ARG 30` a `3.91 Å` e `ASP 65 – LYS 97` a `3.96 Å`).
- **Check della Functional Mask (`3AAA`):** I 24 residui estratti a $< 4.5\text{ Å}$ da CP (`resid 4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36, 41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100`) includono tutti i residui della faccia concava e dei loop $\beta$-hairpin descritti da *Takeda et al. (2010)*.

### ⚠️ Punti critici o di dubbio da evidenziare al Professore
1. **Limite del singolo modello NMR (`2MYO` vs `1MYO`):** `2MYO` è una singola struttura media minimizzata e non rappresenta l'intero ensemble conformazionale NMR (`1MYO` contiene 44 conformeri). Inoltre, l'analisi del `REMARK 500` rivela un contatto ravvicinato tra atomi pesanti del backbone (**`HIS 80 O – SER 82 N` a `2.14 Å`**) e **35 residui fuori dalle regioni di Ramachandran standard** (tra cui **`HIS 68`**, **`ASP 93`** e **`THR 103`**). Questi outlier geometrici da soli non dimostrano instabilità locale, ma richiedono minimizzazione energetica e verifica dinamica in MD.
2. **Limite della Functional Mask computazionale:** Escludere i 24 residui di contatto con CP preserva il sito di legame nel modello, ma solo un test sperimentale *in vitro* può confermare che l'affinità di binding per CP resti inalterata.

---

## 📌 PARTE 2 : Preparazione del Sistema (`CHARMM36m`, Istidine, Solvatazione `TIP3P`, Ionizzazione) e Benchmark

### 🗣️ Parte 2 : Discorso da dire al Professore
> *"Una volta validata la struttura di partenza, nella **Parte 2** abbiamo preparato il sistema Wild-Type per la Dinamica Molecolare in **NAMD 3** con il force field all-atom **CHARMM36m**.*
> *Prima di ricostruire la topologia con `psfgen`, abbiamo assegnato esplicitamente lo stato di protonazione delle **4 Istidine (`HIS 40, 67, 68, 80`)** a pH fisiologico neutro ($\text{pH} \approx 7.0$, dove l'anello imidazolico è neutro ma può esistere come tautomero `HSD` protonato su `ND1` o `HSE` protonato su `NE2`). Ispezionando i donatori e gli accettori entro $3.5\text{ Å}$:*
> - *Abbiamo assegnato **`HIS 40 -> HSD`** perché il suo `ND1-HD1` forma il ponte a idrogeno inter-repeat canonico del motivo `TPLH` con l'ossigeno di `SER 74 OG` ($2.44\text{ Å}$), e abbiamo mantenuto **`HSD`** per le due istidine esposte **`HIS 68`** e **`HIS 80`**;*
> - *Abbiamo invece assegnato **`HIS 67 -> HSE`** (protonata su `NE2`) perché a $3.24\text{ Å}$ dal suo `ND1` si trova il gruppo $-\text{NH}_3^+$ carico positivamente di `LYS 37` (da cui `ND1` deprotonato può accettare un H-bond senza repulsione elettrostatica), mentre `NE2-HE2` dona un H-bond all'ossigeno carbonilico di `GLY 35 O` ($2.40\text{ Å}$).*
> *Abbiamo quindi rimosso gli idrogeni NMR (`protein and noh`), ricostruito gli idrogeni secondo la topologia `CHARMM36m` (ottenendo una carica netta della proteina di **$-5\text{ e}$**), solvatato la proteina in una scatola ortorombica di acqua esplicita **`TIP3P` con padding di $12\text{ Å}$** (`71.84 × 58.97 × 54.84 Å`) e ionizzato a concentrazione fisiologica **$0.15\text{ M NaCl}$** sostituendo 41 molecole d'acqua con **`23` ioni $\text{Na}^+$ e `18` ioni $\text{Cl}^-$** (`6.535` molecole d'acqua, **`21.461` atomi totali**, carica netta exactamente pari a **`0.0000 e`**).*
> *Infine, seguendo la Slide 21 (*Benchmark first*), abbiamo misurato le prestazioni di NAMD 3 sul nostro hardware, trovando che l'uso degli **8 core fisici puri (`namd3 +p8` = `12.68 ns/day`)** è del `+40%` più veloce rispetto a 16 thread logici SMT (`9.01 ns/day`)."*

### ✅ Perché è valido e in che modo oggettivo l'abbiamo checkato
- **Check indipendente della topologia e geometria ([`report/02_scripts/01_validate_wt_preparation.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_scripts/01_validate_wt_preparation.tcl) & [`report/02_utilities/01_wt_preparation_checks.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_utilities/01_wt_preparation_checks.txt)):**
  - Spostamento massimo dei 118 atomi $C_\alpha$ dopo `psfgen`: **`0.0000 Å`** (`900/902` atomi pesanti identici al millesimo di Ångström; solo i 2 ossigeni carbossilici C-terminali `OT1/OT2` di `GLN 118` sono stati riposizionati di `~0.6 Å` dalla patch `CTER` di CHARMM).
  - Verifica neutralità elettrica per il metodo Particle Mesh Ewald (PME): $-5\text{ e (proteina)} + 23\text{ e (Na}^+) - 18\text{ e (Cl}^-) = \mathbf{0.0000\text{ e}}$.
  - Verifica contatti sterici post-solvatazione: **`0` contatti proteina–acqua o proteina–ioni $< 2.0\text{ Å}$**.

### ⚠️ Punti critici o di dubbio da evidenziare al Professore
1. **Tautomeri delle istidine come assunzione modellistica fissa:** L'assegnazione `HSD`/`HSE` si basa sulla geometria locale statica del PDB NMR e non su una titolazione quantistica o su calcoli espliciti di $\text{p}K_a$ ( né può cambiare dinamicamente durante una MD classica a carica fissa).
2. **Margini della box PBC sotto stress termico:** Il padding di $12\text{ Å}$ (con cutoff non-bonded a $12\text{ Å}$) è calcolato sulla conformazione iniziale compatta; se la proteina si srotolasse completamente a 400 K, i margini verso le immagini periodiche diminuirebbero (per questo controlliamo sempre che il raggio di girazione $R_g$ resti contenuto intorno a `14.5–15.7 Å`).

---

## 📌 PARTE 3 : Simulazioni Pilota del Wild-Type (`300 K` vs `400 K`), Corretta Definizione del $Q$-value e Regioni Candidate

### 🗣️ Parte 3 : Discorso da dire al Professore
> *"Nella **Parte 3** abbiamo eseguito le due simulazioni pilota NPT del Wild-Type: **`1.5 ns` a `300 K`** (controllo a temperatura ambiente) e **`2.5 ns` a `400 K`** (stress termico accelerato secondo il protocollo di Bekker et al., 2019), dopo `1.000` step di minimizzazione energetica.*
> *Qui abbiamo affrontato e risolto subito un **punto metodologico importante sul calcolo del $Q$-value (frazione di contatti nativi)**:*
> - *In un primo test avevamo caricato in VMD il file `.pdb` non minimizzato come `Frame 0` prima del `.dcd`. Ci siamo però accorti che durante i `1.000` step di minimizzazione energetica iniziale gli atomi troppo vicini del modello NMR si rilassano leggermente, facendo perdere istantaneamente circa il `17%` dei contatti entro $4.5\text{ Å}$ già prima del primo frame del `.dcd` (per cui nei grafici il $Q$ partiva da `~0.825` invece che da `1.0`, e poi scendeva intorno a `0.70`).*
> - *Allineandoci invece esattamente al Suo script [`analysisMD.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/analysisMD.tcl) (`mol new system.psf` + `mol addfile production.dcd`), abbiamo preso come riferimento nativo (`Frame 0`) il **primo frame della traiettoria `.dcd` (dopo la minimizzazione)**: in questo modo **$Q(t=0)$ parte esattamente da `1.00`** e il **$Q$-value del Wild-Type si mantiene esattamente nella fascia `0.80–0.85`** (e a 300 K intorno a `0.78–0.85`).*
>
> *Confrontando `300 K` e `400 K` attraverso tutte le osservabili (`RMSD`, `Q-value`, `Rgyr`, `SASA`, `% Alpha-Helix`, `H-bonds` e `RMSF`):*
> - *A **`300 K`** la proteina rimane compatta e stabile ($\text{RMSD} \approx 2.6\text{ Å}$, $R_g \approx 14.51\text{ Å}$, $\text{SASA} \approx 7.507\text{ Å}^2$, `14.3` H-bonds medi).*
> - *A **`400 K`** osserviamo un'espansione globale ($\text{RMSD}$ sale a `6.68 Å` nell'ultimo 30%, $R_g$ sale a `15.74 Å`, $\text{SASA}$ sale a `8.840 Å²`), ma il contenuto complessivo di $\alpha$-elica non crolla (`~46–53%`): questo dimostra che a 400 K **le singole eliche non si sciolgono uniformemente, ma cedono i contatti terziari tra i moduli Ankyrin e si aprono i loop di giunzione**!*
> *Dall'analisi combinata di RMSF, perdita di occupazione degli H-bond e ritenzione dei contatti per residuo ($Q_i$ nell'ultimo 30% a 400 K), abbiamo identificato **4 regioni candidate** fuori dalla maschera funzionale:*
> 1. ***Regione 1 (`Residui 28–37`, loop e ingresso del 1° repeat centrale `ANK2`):** qui l'H-bond `LEU 32 – ARG 36` crolla dal `43.6%` (300 K) allo `0.08%` (400 K), e **`LYS 37`** mostra un balzo di RMSF da `0.95 Å` a `2.96 Å` formando un ponte salino distorsivo con `ASP 60`;*
> 2. ***Regione 2 (`Residui 46–55`, core idrofobico del 1° repeat centrale `ANK2`):** i residui sepolti `ILE 50, LEU 51, PHE 53` perdono gran parte dei loro contatti nativi perché l'intero primo repeat centrale è il modulo più vulnerabile della proteina;*
> 3. ***Regione 3 (`Residui 60–69`, giunzione centrale tra 1° e 2° repeat, `ASP 60 – HIS 68`):** qui crollano gli H-bond `THR 70 – ILE 61` (`24.6% -> 0.08%`) e `LYS 97 – ASP 65` (`47.9% -> 5.4%`), con forti fluttuazioni su `ASN 62 – PRO 64` e tensione sul $\beta$-turn di **`HIS 68`**;*
> 4. ***Regione 4 (`Residui 94–103`, ingresso della regione C-terminale `ANK4`):** qui si perde il legame a idrogeno di capping **`THR 103 OG1 – ASP 100 O`** (`15.56% -> 0.32%`), innescando lo sfilacciamento C-terminale (`LEU 102 – THR 103` con $\text{RMSF} > 3.1\text{ Å}$)."*

### ✅ Perché è valido e in che modo oggettivo l'abbiamo checkato
- **Cross-check tra 7 osservabili indipendenti ([`scripts/03_analyze_wt_baseline.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/03_analyze_wt_baseline.tcl) & [`report/02_wild_type_characterization_and_baseline.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/report/02_wild_type_characterization_and_baseline.md)):** Non abbiamo scelto i siti guardando solo l'RMSF (come avverte la Slide 22: *"flexibility alone does not establish local instability"*), ma abbiamo verificato che nelle 4 regioni candidate coincidessero **alta fluttuazione RMSF + basso $Q_i$ locale + rottura misurata di specifici legami a idrogeno (`hbonds_details_300K.txt` vs `400K.txt`)**.
- **Scomposizione del $Q$-value su 3 livelli:** Globale/Per classe chimica di Bekker (`all-all`, `hydrophilic-all`, `hydrophobic-small`), per modulo strutturale (`ANK1..4`) e per singolo residuo ($Q_i$).

### ⚠️ Punti critici o di dubbio da evidenziare al Professore
1. **Assenza di uno stadio separato di pre-equilibrazione nelle run pilota:** Nelle simulazioni pilota la dinamica NPT è partita direttamente dopo i `1.000` step di minimizzazione (senza uno step preliminare con *position restraints* sul backbone a 300 K). Per ridurre l'effetto del rilassamento iniziale abbiamo calcolato i $Q_i$ solo sull'**ultimo 30% della traiettoria (`1.75–2.50 ns`)**, ma parte della discesa iniziale del $Q$ nei primi picosecondi (comune a 300 K e 400 K) riflette l'adattamento termico iniziale.
2. **Durate diverse (`1.5 ns` a 300 K vs `2.5 ns` a 400 K) e singola replica esplorativa:** Per confrontare in modo equo 300 K e 400 K nel report abbiamo confrontato sia la finestra comune (`0–1.5 ns`) sia l'ultimo 30% (`1.75–2.50 ns`), trattando queste run come **traiettorie esplorative** per guidare lo screening mutazionale.

---

## 📌 PARTE 4 : Screening Mutazionale Multi-Piattaforma e Modellazione 3D (`FoldX BuildModel`)

### 🗣️ Parte 4 : Discorso da dire al Professore
> *"Individuate le 4 regioni candidate nella Parte 3, nella **Parte 4** abbiamo progettato e filtrato le mutazioni stabilizzanti.*
> *Escludendo i 24 residui della Functional Mask (`3AAA`) e i 6 residui delle code terminali (`1–3` e `116–118`), abbiamo sottoposto tutte le **88 posizioni modificabili dello scaffold ($88 \times 20 = 1.760$ varianti)** a uno screening di consenso basato su **4 approcci ortogonali**:*
> 1. ***FoldX 5** (seguendo il workflow della Slide 19: `RepairPDB` che ha ridotto l'energia di clash da `69.1` a `13.5 kcal/mol` $\rightarrow$ `PositionScan` completo sulle 88 posizioni $\rightarrow$ `BuildModel` con 5 run indipendenti e *matched WT references*);*
> 2. ***FireProt 2.0** (che combina calcoli energetici FoldX/Rosetta con l'analisi evolutiva *Back-to-Consensus* su Multiple Sequence Alignment);*
> 3. ***ThermoMPNN** (Graph Neural Network addestrata su dataset sperimentali di $T_m$ per mutazioni singole e doppie);*
> 4. *L'**allineamento al consenso universale degli Ankyrin Repeats** di Mosavi et al. (2002).*
>
> *Incrociando i 4 metodi con i punti deboli della MD a 400 K, abbiamo scartato i falsi positivi (come `D60L`, che aveva un buon $\Delta\Delta G$ in FoldX ma introduceva una Leucina idrofobica sulla superficie esposta di un loop polare ed era infatti bocciato da FireProt) e abbiamo selezionato **3 mutazioni chiave**, tutte e tre di tipo **Back-to-Consensus** e ciascuna mirata a una delle regioni critiche della Parte 3:*
> 1. ***`K37L` (nella Regione 1/2, 1° repeat centrale `ANK2` — 5/5 pilastri):** sostituisce `LYS 37` con la **Leucina canonica in Posizione 4 del motivo consenso `TPLH`** (dove `ANK1`, `ANK3` e `ANK4` hanno già una Leucina!). Elimina il ponte salino distorsivo con `ASP 60` e compatta il core idrofobico di `ANK2`. In `FoldX BuildModel` ($n=5$) dà **$\Delta\Delta G = -2.44 \pm 0.87\text{ kcal/mol}$** (confermato da FireProt `-1.61 kcal/mol` e ThermoMPNN `-0.84 kcal/mol`).*
> 2. ***`T103L` (nella Regione 4, ingresso C-cap `ANK4` — 5/5 pilastri):** ripristina anche qui la **Leucina canonica del consenso Ankyrin** al posto di `THR 103`, rinforzando il packing idrofobico C-terminale ($\Delta E_{\text{VdW+SolvH}} = -4.03\text{ kcal/mol}$, $\Delta\Delta G_{\text{BuildModel}} = -1.15 \pm 0.14\text{ kcal/mol}$, FireProt `-1.78 kcal/mol`, ThermoMPNN `-0.55 kcal/mol`).*
> 3. ***`H68G` (nella Regione 3, giunzione centrale `ANK2–ANK3` — 4/5 pilastri):** è la **#1 scelta Back-to-Consensus di FireProt 2.0** (la Glicina è conservata nel **78%** delle sequenze omologhe) e la **#1 assoluta di ThermoMPNN** (`-1.44 kcal/mol`), con $\Delta\Delta G_{\text{BuildModel}} = -1.73 \pm 0.01\text{ kcal/mol}$. Sostituisce `HIS 68` (outlier di Ramachandran sul $\beta$-turn centrale) con una Glicina che rilassa la tensione sterica di cerniera.*
>
> *Costruendo i mutanti combinati in `FoldX BuildModel` ($n=5$), abbiamo verificato che il **Doppio Mutante `K37L + T103L`** è perfettamente additivo (**$\Delta\Delta G = -3.65 \pm 0.66\text{ kcal/mol}$**, $\Delta\Delta\Delta G_{\text{coupling}} = -0.06\text{ kcal/mol}$) e che il **Triplo Mutante `K37L + H68G + T103L`** raggiunge **$\Delta\Delta G = -5.03 \pm 0.21\text{ kcal/mol}$** (`Best = -5.39 kcal/mol`)."*

### ✅ Perché è valido e in che modo oggettivo l'abbiamo checkato
- **Consenso multi-algoritmo ([`03_mutations/consensus_master_table.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/consensus_master_table.txt) & [`03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/ANALISI_MUTAZIONI_FIREPROT_E_CONSENSO.md)):** Richiedere la concordanza simultanea di FoldX (empirico), Rosetta/FireProt (fisico + MSA), ThermoMPNN (Deep Learning) e Consenso Ankyrin abbatte drasticamente il tasso di falsi positivi dei singoli predittori (`Lecture05` Slide 57–68).
- **Verifica di additività/epistasi con *Matched WT References* ([`03_mutations/buildmodel_epistasis_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/03_mutations/buildmodel_epistasis_report.txt)):** In `FoldX BuildModel` ogni mutante singolo, doppio e triplo è stato valutato su 5 run indipendenti sottraendo l'energia del rispettivo `WT_2MYO_Repair_*.pdb` ottimizzato nello stesso intorno locale (Slide 19).

### ⚠️ Punti critici o di dubbio da evidenziare al Professore
1. **Approssimazione a backbone fisso e solvente implicito di FoldX (`Slide 17`):** Come sottolineato nella Slide 17 (*"FoldX predicts an energetic effect. MD tests structural behaviour under chosen conditions"*), FoldX mantiene il backbone quasi rigido: non può dirci se la rottura del ponte salino `LYS 37 – ASP 60` in `K37L` renda più flessibile il loop `60–68` in acqua calda. Solo la MD esplicita a 400 K può verificarlo.

---

## 📌 PARTE 5 : Screening MD Comparativo a `400 K` (`WT` vs `K37L`, `T103L`, `Doppio`, `Triplo`) e Scoperta Cooperativa

### 🗣️ Parte 5 : Discorso da dire al Professore
> *"Nella **Parte 5** abbiamo sottoposto i 4 candidati (`K37L`, `T103L`, `Doppio K37L+T103L` e `Triplo K37L+H68G+T103L`) a un primo ciclo di **Matched MD Comparison a `400 K` per `2.5 ns` (`Seed 1`)**, preparando tutti i sistemi con lo stesso identico protocollo del Wild-Type (`CHARMM36m + TIP3P 12 Å + 0.15 M NaCl`, stessa box PBC `72.59 × 58.52 × 58.88 Å`, carica totale `+0.0000 e`).*
> *I risultati dinamici a 400 K hanno confermato la stabilizzazione prevista e ci hanno permesso di capire un **meccanismo cooperativo fondamentale tra le 3 mutazioni**:*
> 1. ***Effetto di `K37L` (Singolo 1) e `T103L` (Singolo 2):** Entrambi migliorano nettamente la ritenzione dei contatti nativi rispetto al Wild-Type (`Q: all-all` e `Q: hydrophilic-all` salgono di `+2% / +4%`), riducono la SASA di oltre **`-870 Å²`** e abbassano il Backbone RMSD da `6.68 Å` (WT) a **`4.17 Å` (`K37L`)** e **`3.77 Å` (`T103L`)**. In particolare, in `K37L` i residui del core idrofobico sepolto di `ANK2` (`ILE 50, LEU 51, PHE 53`, che avevamo identificato nella Parte 3) **recuperano tutti i loro contatti nativi**!*
> 2. ***Il piccolo effetto collaterale sul loop `ASP 60` in `K37L` e nel Doppio `K37L+T103L`:** Nel Doppio `K37L+T103L` raggiungiamo la massima compattezza globulare ($R_g = 14.76\text{ Å}$ contro `15.74 Å` del WT) e un ottimo $Q$-value; tuttavia, eliminando `LYS 37`, il suo ex-partner **`ASP 60`** (sulla punta del loop centrale `60–68` tra `ANK2` e `ANK3`) rimane privo di ancoraggio e oscilla maggiormente nel solvente ($\text{RMSF}_{60}$ sale da `2.53 Å` a `5.00 Å` in `K37L` e `3.35 Å` nel Doppio), mantenendo alto l'RMSD globale del Doppio (`6.40 Å`).*
> 3. ***Perché il Triplo Mutante (`K37L + H68G + T103L`) risolve il problema e vince su tutte le metriche:** Aggiungendo **`H68G`** esattamente al residuo `68` (l'uscita del loop `60–68` tra i due repeat centrali), la Glicina 68 rilassa la tensione del $\beta$-turn e **blocca la fluttuazione del loop `ASP 60`, facendo scendere il suo RMSF da `5.00 / 3.35 Å` a soli `1.50 Å`** (e quello di `LEU 37` a **`0.91 Å`** a 400 K!). Di conseguenza, il **Triplo Mutante** ottiene contemporaneamente:*
>    - *Il **miglior $Q$-value globale** e sul modulo critico `ANK2`;*
>    - *Il **minimo Backbone RMSD globale (`3.30 ± 0.42 Å` contro `6.68 Å` del WT, `-51%`)**;*
>    - *La **minima deformazione del sito di legame all'actina (`3AAA` RMSD = `2.83 ± 0.26 Å` contro `7.27 Å` del WT, `-61%`)**;*
>    - *La **minima SASA totale (`7.647 Å²`, `-1.192 Å²` rispetto al WT)** e il minor numero assoluto di residui instabili (`-44%` rispetto al WT)."*

### ✅ Perché è valido e in che modo oggettivo l'abbiamo checkato
- **Check di qualità topologica Gate 2 ([`04_md_simulations/matched_systems_qc_report.txt`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/04_md_simulations/matched_systems_qc_report.txt)):** Verificato che tutte le varianti con `K37L` (carica proteina `-6 e`) abbiano esattamente `26 Na+ / 20 Cl-` e quelle neutre (`-5 e`) abbiano `25 Na+ / 20 Cl-`, garantendo $q_{\text{tot}} = +0.0000\text{ e}$ e volumi di cella PBC identici al millesimo di Ångström (`72.589 × 58.517 × 58.883 Å`).
- **Analisi sistematica separata e unificata ([`05_analysis/single_runs/`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/05_analysis/single_runs) & [`06_presentation/fig9..11.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/06_presentation)):** Ogni traiettoria ha `1.250` frame completi (`2.50 ns`) ed è stata valutata nella stessa finestra finale del 30% (`1.75–2.50 ns`).

### ⚠️ Punti critici o di dubbio da evidenziare al Professore (che motivano la Parte 6!)
1. **Singola replica (`Seed 1`) nello screening a 5 sistemi:** Avendo testato 5 sistemi diversi per confrontare singoli, doppio e triplo, in questa fase di screening abbiamo eseguito 1 solo velocity seed (`12345`) per sistema. Poiché la dinamica molecolare a 400 K è stocastica e caotica, una singola run non permette di stimare l'incertezza statistica ($\mu \pm \sigma$) tra repliche indipendenti.
2. **Assenza di pre-equilibrazione restrained a 300 K prima del run a 400 K:** Per mantenere il protocollo 100% *matched* con la baseline WT della Parte 3, anche le 4 run mutanti sono partite direttamente a 400 K dopo i 1.000 step di minimizzazione.

---

## 📌 PARTE 6 : Proposta della Fase Finale di Validazione da Condividere con il Professore

### 🗣️ Parte 6 : Discorso da dire al Professore
> *"Professore, alla luce di questi risultati intermedi e delle due osservazioni metodologiche che abbiamo appena evidenziato (singolo velocity seed e assenza di uno step separato di pre-equilibrazione a 300 K prima dello stress a 400 K), per chiudere il progetto in modo rigoroso prima della presentazione del 16 Dicembre vorremmo procedere così:*
> 1. ***Concentrare la validazione finale sul confronto diretto tra il `Wild-Type` e il nostro candidato vincitore (`Triplo Mutante K37L+H68G+T103L`, ed eventualmente il `Doppio K37L+T103L`)**, dato che lo screening sul Seed 1 ha già chiarito il ruolo meccanico dei singoli mutanti.*
> 2. ***Introdurre nel protocollo NAMD uno stadio esplicito di Pre-Equilibrazione a `300 K` (`100 ps` in NPT a `300 K` con vincoli armonici sul backbone, `constraints on`)** subito dopo i `1.000` step di minimizzazione, in modo da far rilassare perfettamente l'acqua `TIP3P`, gli ioni e le catene laterali a temperatura fisiologica prima di rilasciare i vincoli e avviare la produzione a **`400 K` (`2.5 ns`)**.*
> 3. ***Eseguire $n = 3$ repliche indipendenti con `velocity seeds` diversi** per il Wild-Type e per il mutante vincitore con questo protocollo pre-equilibrato, così da presentare nella tabella finale le **medie $\pm$ deviazioni standard ($\mu \pm \sigma$)** e verificare che il guadagno di $Q$-value e la riduzione dell'RMSD siano statisticamente significativi e indipendenti dal seme casuale.*
>
> *Secondo Lei questa strategia per la validazione finale è corretta o ci consiglia di modificare/aggiungere qualche controllo?"*

---

## 🗺️ Tabella Riassuntiva Rapida ("Cheat Sheet" da tenere sotto mano durante il colloquio)

| Fase / Parte | Cosa abbiamo fatto | Risultato / Scoperta Chiave | Come lo abbiamo checkato | Punto Critico & Soluzione |
| :--- | :--- | :--- | :--- | :--- |
| **Parte 1**<br>*(Struttura & Funzione)* | Analisi `2MYO` (118 aa, NMR) + Allineamento su `3AAA` (complesso con CP) | `0` atomi mancanti, `3` Cys ridotte ($>9.4\text{ Å}$), **Functional Mask di 24 residui** ($<4.5\text{ Å}$ da CP) | `01_check_missing_residues.tcl`, `01_inspect_2myo.tcl`, SASA/RSA (`19` res. sepolti) | `2MYO` è 1 solo modello NMR con 35 outlier Ramachandran (`HIS68, THR103`) $\rightarrow$ serve MD |
| **Parte 2**<br>*(Setup CHARMM36m)* | Tautomeri His (`40,68,80=HSD`, `67=HSE`), `TIP3P` $12\text{ Å}$, $0.15\text{ M NaCl}$ | `21.461` atomi, carica netta `0.0000 e`; Benchmark `+p8` = **`12.68 ns/day`** (`+40%` vs `+p16`) | `01_validate_wt_preparation.tcl` ($\Delta C_\alpha = 0.000\text{ Å}$, `0` clash $<2.0\text{ Å}$) | Istidine assegnate da geometria H-bond senza $\text{p}K_a$ esplicito |
| **Parte 3**<br>*(Baseline WT 300K/400K)* | MD pilota `300 K` (`1.5 ns`) vs `400 K` (`2.5 ns`) + correzione riferimento $Q(0)=1.0$ | A 400 K la proteina si espande ($R_g = 15.74\text{ Å}, \text{RMSD} = 6.68\text{ Å}$). Individuate **4 regioni deboli**: `K37`, core `I50-L51`, giunzione `D60-H68`, C-cap `L102-T103` | Riferimento `Frame 0` = 1° frame DCD (`analysisMD.tcl`) $\rightarrow$ $Q(0)=1.0$ e $Q_{\text{WT}} \in [0.80, 0.85]$ | Nessuna pre-equilibrazione separata a 300 K nelle run pilota $\rightarrow$ risolta nella Parte 6 |
| **Parte 4**<br>*(Design & BuildModel)* | Screening `88` posizioni (`1.760` mut) con `FoldX + FireProt + ThermoMPNN + Consenso` | 3 mutazioni Back-to-Consensus: **`K37L`** (`-2.44`), **`T103L`** (`-1.15`), **`H68G`** (`-1.73`); **Doppio** `-3.65`, **Triplo** **`-5.03 kcal/mol`** | Consenso a 5 pilastri + `FoldX BuildModel` ($n=5$ run con *matched WT references*) | FoldX è a backbone quasi rigido e solvente implicito $\rightarrow$ serve verifica MD a 400 K |
| **Parte 5**<br>*(Screening MD 400K)* | MD `400 K` (`2.5 ns`, `Seed 1`) su `WT`, `K37L`, `T103L`, `Doppio`, `Triplo` | `K37L` salva `ANK2` ma lascia mobile `ASP 60`; aggiungendo **`H68G` nel Triplo**, $\text{RMSF}_{60}$ scende a **`1.50 Å`**, $\text{RMSD} = \mathbf{3.30\text{ Å}}$, sito `3AAA` = **`2.83 Å`**, `-44%` res. instabili! | Gate 2 QC (`matched_systems_qc_report.txt`) + 15 plot singoli e 3 plot unificati (`fig9–11`) | Screening eseguito su 1 solo seed senza pre-equilibrazione a 300 K $\rightarrow$ porta alla Parte 6 |
| **Parte 6**<br>*(Validazione Finale)* | Proposta Passo B su `WT` vs `Vincitore (Triplo / Doppio)` | Pre-equilibrazione `100 ps` a `300 K` (`constraints on` sul backbone) + `2.5 ns` a `400 K` su **$n=3$ velocity seeds** ($\mu \pm \sigma$) | Protocollo multi-stadio `Lecture04` Slide 6–15 + test statistico $t$ di Welch | Chiude tutti i punti aperti per la presentazione del 16 Dicembre |
