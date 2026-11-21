# 🧬 PUNTO 1 DELLA SCALETTA DEL PROFESSORE — "1. Analyze the assigned protein: Function, chain and biological assembly. Missing residues or atoms"

**Riferimento Ufficiale nella Scaletta del Corso:**  
- [`Laboratory_Project_2026.pdf`](file:///home/olly/UNI/BioMolEng2026/laboratorio/Laboratory_Project_2026.pdf) — **Slide 3** (*Learning Objective 1: Protein analysis — Connect sequence and structure to the biological function of the assigned protein*), **Slide 4** (*Project outputs: Introduction on the target protein*), **Slide 7** (*Protein Engineering: improve the protein stability without affecting its function*), **Slide 14** (*Protein assignments: Team #10 Nedina Oleg & Rigoni Alessia $\rightarrow$ Myotrophin, PDB ID `2MYO`, Length 118*), **Slide 15** (*Project workflow — Step 1: Analyze the assigned protein. Function, chain and biological assembly. Missing residues or atoms*), **Slide 19** (*Selection for MD: Combine the energetic prediction with structural plausibility and functional constraints*).

---

## 1. 📖 Introduzione al Problema e Riferimenti alle Slide di Teoria del Corso

In qualsiasi progetto di modellazione biomolecolare e ingegneria proteica (*Protein Engineering*), il primo passo non è mai lanciare un algoritmo alla cieca, ma comprendere a fondo la **chimica, la struttura tridimensionale e la funzione biologica** della proteina assegnata, verificando l'integrità del file PDB di partenza (*"Garbage In, Garbage Out"*):

1. **Legame tra Struttura Primaria, Secondaria, Terziaria e Funzione (`Lecture01 - Introduction.pdf`, Slide 13–48):**
   - Le proteine sono polimeri lineari dei 20 amminoacidi naturali uniti da legami peptidici planari (`Lecture01`, Slide 14–26: *Primary structure & Peptidic bond*), divisi per proprietà chimico-fisiche in **idrofobici** (che formano il core interno lontano dall'acqua, Slide 17–18), **polari** (che formano legami a idrogeno, Slide 19) e **carichi** ($\pm$, che formano ponti salini a pH neutro, Slide 20–21).
   - La catena si ripiega in elementi di **struttura secondaria** ($\alpha$-eliche стабилиizzate da legami a idrogeno backbone $i \to i+4$ con passo di $3.6$ residui/giro, Slide 27–30, e $\beta$-sheet/turn, Slide 31–34) che si impacchettano nella **struttura terziaria** globulare o modulare (Slide 36–38) ed eventualmente in complessi **quaternari** (Slide 39).
   - **Il vincolo fondamentale dell'ingegneria proteica (`Laboratory_Project_2026.pdf`, Slide 7 & Slide 19; `Lecture05`, Slide 58):** quando progettiamo mutazioni per innalzare la temperatura di melting ($T_m$, `Laboratory_Project_2026.pdf` Slide 5–6), **non dobbiamo mai mutare i residui responsabili della funzione biologica** (sito catalitico o interfaccia di legame proteina–proteina).
2. **Origine Sperimentale del File PDB e Preparazione della Configurazione Iniziale (`Lecture02 - Protein Structure Determination.pdf`, Slide 2–3, 25–43; `Lecture03 - Molecular Dynamics Setup.pdf`, Slide 2, 19–25):**
   - Prima di impostare una simulazione di Dinamica Molecolare, occorre selezionare e ispezionare la configurazione atomica iniziale scaricata dal **Protein Data Bank (PDB)** (`Lecture02`, Slide 3, 40–42; `Lecture03`, Slide 2).
   - È essenziale distinguere il metodo sperimentale (`EXPDTA` nel PDB, `Lecture02` Slide 25–26):
     - Le strutture da **Cristallografia a Raggi X** spesso presentano loop flessibili non risolti (*missing residues / gaps* nella catena), atomi di catene laterali mancanti, molecole d'acqua o co-soluti di cristallizzazione (`HOH`, `HETATM`) e copie simmetriche di impaccamento cristallino, ma **non contengono atomi di idrogeno** (perché l'idrogeno ha un solo elettrone e non diffrange i raggi X, `Lecture02` Slide 27–38).
     - Le strutture da **Risonanza Magnetica Nucleare in Soluzione (Solution NMR)**, invece, vengono risolte in acqua a temperatura e pH noti (`REMARK 210`), hanno la catena polipeptidica completa e **contengono già gli atomi di idrogeno sperimentali**, ma con una nomenclatura NMR (`H1, H2, HB2, HB3...`) incompatibile con il generatore topologico `psfgen` di CHARMM, oltre a presentare talvolta contatti ravvicinati (*close contacts*) o angoli diedri $(\phi, \psi)$ di Ramachandran leggermente tesi (`REMARK 500`).
3. **Ponti Disolfuro e Stati di Protonazione a pH Fisiologico (`Lecture01`, Slide 21 & 38; `Lecture03`, Slide 2 & 15):**
   - Due cisteine (`CYS`) spazialmente vicine ($d(S_\gamma\text{--}S_\gamma) \approx 2.05\text{ Å}$) possono ossidarsi formando un legame covalente disolfuro (`-S-S-`, patch `DISU` in CHARMM, `Lecture01` Slide 38).
   - L'istidina (`HIS`) ha un $pK_a \approx 6.0\text{--}6.5$ vicino al pH fisiologico: nel force field CHARMM non esiste un residuo generico `HIS`, ma bisogna specificare esplicitamente quale dei tre stati chimici adottare (`Lecture03`, Slide 15):
     - **`HSD`**: istidina neutra con il protone sull'azoto $N_{\delta 1}$ (`ND1-HD1`) e doppietto elettronico libero (accettore H-bond) su $N_{\epsilon 2}$ (`NE2`);
     - **`HSE`**: istidina neutra con il protone sull'azoto $N_{\epsilon 2}$ (`NE2-HE2`) e doppietto libero su $N_{\delta 1}$ (`ND1`);
     - **`HSP`**: istidina doppiamente protonata carica $+1$ (entrambi gli azoti protonati, tipica di $\text{pH} < 6.0$).

---

## 2. 🎯 Cosa Cercavamo (Le Domande Scientifiche del Punto 1)

Per soddisfare con rigore da ricerca il **Punto 1** della scaletta del Prof. Gautieri sulla nostra proteina assegnata **Myotrophin (`2MYO`)**, dovevamo rispondere senza alcuna ipotesi a 5 domande precise:
1. **Identità e Funzione Biologica:** Qual è il ruolo fisiologico della Myotrophin (nota anche come proteina **V-1**), da quale organismo proviene `2MYO` e qual è la sua architettura strutturale (fold, domini ripetuti, struttura secondaria)?
2. **Qualità e Completezza del File PDB (`2MYO.pdb`):** Qual è il metodo sperimentale? Quanti modelli e quante catene sono presenti? L'unità biologica attiva (*biological assembly*) è un monomero o un oligomero? Ci sono residui o atomi pesanti mancanti nella catena `1..118`? Ci sono discrepanze di numerazione rispetto a UniProt?
3. **Stato Chimico delle Cisteine (`CYS`):** Quante cisteine ci sono in `2MYO` e a che distanza reciproca si trovano i loro atomi di zolfo $S_\gamma$? Formano ponti disolfuro o sono tioli liberi ridotti (`-SH`)?
4. **Assegnazione Razionale dei Tautomeri delle Istidine (`HIS`):** Quante istidine ci sono in `2MYO`? Quali donatori e accettori di legame a idrogeno si trovano entro $3.5\text{ Å}$ dai loro anelli imidazolici (`ND1` e `NE2`)? Quale tautomero (`HSD` vs `HSE` vs `HSP`) preserva la rete nativa di legami a idrogeno evitando repulsioni sterico-elettrostatiche?
5. **Mappatura Oggettiva del Sito Funzionale (*Functional Masking*):** Quali amminoacidi della Myotrophin entrano in contatto diretto con il suo bersaglio biologico (**Actin-Capping Protein, CP**) e devono quindi essere **tassativamente esclusi (mascherati)** da qualsiasi mutazione nella Fase 3?

---

## 3. 🛠️ Come L'Abbiamo Cercato (Metodologia, File e Script Eseguiti)

Per rispondere alle 5 domande senza mai tirare a indovinare ([`AGENTS.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/AGENTS.md) Regola 3: *Zero Ipotesi*), abbiamo combinato l'analisi documentale dei record PDB con misure geometriche quantitative in **VMD 2.0.0** ([`scripts/01_inspect_2myo.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/01_inspect_2myo.tcl)):

1. **Acquisizione delle Strutture di Riferimento (`01_structures/raw/`):**
   - Scaricato [`01_structures/raw/2MYO.pdb`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/2MYO.pdb) (struttura NMR in soluzione della Myotrophin di *Rattus norvegicus*, *Yang et al., Structure 1998*).
   - Scaricato [`01_structures/raw/3AAA_complex_ref.pdb`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/01_structures/raw/3AAA_complex_ref.pdb) (struttura cristallografica a Raggi X a risoluzione $2.20\text{ Å}$ del complesso biologico tra l'eterodimero $\alpha/\beta$ dell'**Actin-Capping Protein**, `Chain A` e `Chain B`, e la **Myotrophin/V-1**, `Chain C`, *Takeda et al., PLoS Biology 2010*).
2. **Ispezione dei Record Ufficiali di `2MYO.pdb` (righe `1–200`):**
   - Letti direttamente i campi `HEADER`, `TITLE`, `COMPND`, `SOURCE`, `EXPDTA`, `REMARK 210` (condizioni sperimentali NMR), `REMARK 215`, `REMARK 300/350` (stato oligomerico `BIOMT`), `REMARK 500` (*close contacts* e outlier torsionali di Ramachandran), `DBREF` (allineamento con UniProt `P62775`), `SEQRES` e `HELIX 1–8`.
3. **Esecuzione dello Script Automatico di Ispezione in VMD ([`scripts/01_inspect_2myo.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/01_inspect_2myo.tcl)):**
   - **Verifica completezza e FASTA:** conteggio di `[atomselect top "protein"]`, `[atomselect top "protein and noh"]`, `[atomselect top "protein and name CA"]` ed estrazione della sequenza continua tramite `getFasta` (da [`scripts/VMDextensions.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/VMDextensions.tcl)).
   - **Misura delle distanze $S_\gamma\text{--}S_\gamma$ delle Cisteine:** selezione `resname CYS and name SG` e calcolo di `vecdist` su tutte le coppie $(i, j)$.
   - **Scansione della sfera di coordinazione delle 4 Istidine:** per ciascuna istidina (`resid 40, 67, 68, 80`) e per ciascuno dei due azoti dell'anello (`ND1` e `NE2`), selezione automatica `(noh and not resid $hres) and within 3.5 of (resid $hres and name $ring_n)` e calcolo delle distanze esatte al centesimo di Ångström.
   - **Sovrapposizione Strutturale `2MYO` vs `3AAA` ed Estrazione della Functional Mask:**
     - Allineamento roto-traslazionale (`measure fit`) dei carboni $\alpha$ (`resid 2 to 118 and name CA`) di `2MYO` (`Chain A`) sulla Myotrophin del complesso `3AAA` (`Chain C`) e calcolo di `measure rmsd`.
     - Estrazione oggettiva di tutti i residui di Myotrophin aventi almeno un atomo pesante a distanza **$< 4.5\text{ Å}$** dalle subunità $\alpha$ e $\beta$ dell'Actin-Capping Protein:
       ```tcl
       set interface_sel [atomselect 1 "protein and chain C and noh and within 4.5 of (protein and chain A B and noh)"]
       ```
4. **Allineamento Strutturale dei 4 Ankyrin Repeats contro il Consenso Universale ([`scripts/04_ankyrin_consensus_and_targets.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/04_ankyrin_consensus_and_targets.py)):**
   - Mappatura posizione per posizione (`1..33`) dei 4 moduli `ANK1–ANK4` contro la firma consenso di *Mosavi, Minor & Peng (PNAS 2002)* (`GXTPLHLAAXXGHLEIVEVLLKXGADVNAXDXX`).

---

## 4. 💡 Cosa Abbiamo Trovato e Considerazioni Critiche (Risultati Quantitativi)

Ecco tutte le evidenze numeriche estratte e documentate in [`STATUS/discovery.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/STATUS/discovery.md) e [`ROADMAP/ROAD1.md`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/ROADMAP/ROAD1.md):

### 4.1 Identikit Strutturale e Completezza al 100% di `2MYO.pdb`
- **Metodo e Condizioni:** Struttura **Solution NMR** (*minimized average structure* derivata da 90 conformeri, `REMARK 210` righe 48–67), acquisita in $\text{H}_2\text{O}$ a **$T = 297\text{ K}$ ($24^\circ\text{C}$)**, **$\text{pH} = 6.2$**, tampone fosfato $50\text{ mM}$.
- **Stato Biologico (`REMARK 350`, riga 93):** **Monomero** (`AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC`, singola catena `Chain A`).
- **Completezza della Catena (`0` gap):**
  - Numero di residui ($C_\alpha$): **`118`** continui da `MET 1` a `GLN 118` (`getFasta` restituisce l'intera stringa `MCDKEFMWALKNGDLDEVKDYVAKGEDVNRTLEGGRKPLHYAADCGQLEILEFLLLKGADINAPDKHHITPLLSAVYEGHVSCVKLLLSKGADKTVKGPDGLTALEATDNQAIKALLQ` senza alcuna interruzione `/`).
  - Numero di atomi nel PDB: **`1.815` atomi totali**, composti da **`902` atomi pesanti (`noh`)** (incluso l'ossigeno C-terminale `OXT`) e **`913` atomi di idrogeno NMR**. Zero molecole d'acqua cristallografica e zero ligandi estranei (`HETATM`).
  - **Nota critica sulla numerazione (`DBREF`, riga 163):** `DBREF 2MYO A 2 118 UNP P62775 MTPN_RAT 1 117`. Nella proteina matura *in vivo* la Metionina iniziale (`MET 1`) viene scissa, quindi `CYS 2` del PDB corrisponde al residuo `1` di UniProt. Averlo verificato subito ci ha evitato qualsiasi errore di "shift di +1" nel leggere gli articoli di letteratura!

### 4.2 Architettura a 4 Ankyrin Repeats (`ANK1–ANK4`) e Confronto Ratto (`2MYO`) vs Uomo (`3AAA`)
- La Myotrophin ha un fold a **solenoide $\alpha$-elica leggermente curvo** (~$46\text{ Å} \times 21\text{ Å}$) formato dall'impilamento lineare di **4 Ankyrin-like Repeats**:
  1. **`ANK1` (N-cap, residui `1–34`):** Elica interna `HELIX 1` (`ASP 3 – LYS 11`), Elica esterna `HELIX 2` (`LEU 15 – VAL 22`), loop inter-repeat `23–34`.
  2. **`ANK2` (Core centrale, residui `35–67`):** Motivo `TPLH`-like (`37–40`), Elica interna `HELIX 3` (`HIS 40 – CYS 45`), Elica esterna `HELIX 4` (`ILE 50 – LEU 56`), loop `GADVNA` con breve elica `HELIX 5` (`ASP 60 – ALA 63`).
  3. **`ANK3` (Core centrale, residui `68–100`):** $\beta$-turn (`68–70`), Elica interna `HELIX 6` (`PRO 71 – GLU 78`), Elica esterna `HELIX 7` (`SER 82 – LEU 88`), loop `GADKTV` (`91–100`).
  4. **`ANK4` (C-cap, residui `101–118`):** Elica C-cap (`LEU 102 – ASN 110`) e coda C-terminale con elica $3_{10}$ `HELIX 8` (`LYS 114 – LEU 116`).
- **Scoperta evolutiva straordinaria (`2MYO` Ratto vs `3AAA` Uomo):**
  - Confrontando `2MYO` (Ratto, NMR libera in soluzione) con la catena `C` di `3AAA` (Uomo, Raggi X legata all'Actin-Capping Protein), **117 amminoacidi su 118 sono identici ($99.2\%$ di identità)**! L'unica differenza in tutta la proteina è al residuo `105` nel C-cap: **`LEU 105` nel Ratto vs `PHE 105` nell'Uomo**.
  - L'**$\text{RMSD}_{C_\alpha}$ (`resid 2–118`) tra la forma libera NMR (`2MYO`) e la forma legata X-Ray (`3AAA`) è di appena `2.365 Å`**, dimostrando che la Myotrophin funziona come un "corpo rigido pre-organizzato" e validando l'eccellente qualità geometrica della nostra struttura NMR di partenza!

### 4.3 Verifica Oggettiva delle 3 Cisteine e Scoperta sulla Rete H-Bond delle 4 Istidine
1. **Le 3 Cisteine (`CYS 2`, `CYS 45`, `CYS 83`) non formano ponti disolfuro:**
   - Distanze misurate: $d(\text{CYS } 2, \text{CYS } 45) = \mathbf{21.88\text{ Å}}$, $d(\text{CYS } 2, \text{CYS } 83) = \mathbf{27.01\text{ Å}}$, $d(\text{CYS } 45, \text{CYS } 83) = \mathbf{9.43\text{ Å}}$ (molto superiori ai $\sim 2.05\text{ Å}$ di un legame disolfuro). Tutte e 3 sono cisteine ridotte (`-SH`), coerente con l'ambiente intracellulare riducente del citoplasma cardiaco.
2. **La Rete di Legami a Idrogeno di `HIS 40` (`HSD`) e `HIS 67` (`HSE`):**
   - Nei tutorial standard si applica spesso `pdbalias residue HIS HSE` a tutte le istidine indistintamente. La nostra ispezione a $3.5\text{ Å}$ ha rivelato che farlo su `2MYO` avrebbe introdotto due errori chimici gravi:
     - **`HIS 40` (nel motivo firma `TPLH` di `ANK2`) $\rightarrow$ `HSD`:** il suo azoto **`ND1`** si trova a soli **`2.44 Å`** dall'ossigeno idrossilico **`OG` di `SER 74`** (nell'elica interna di `ANK3`), mentre `NE2` punta verso `ARG 36` ($3.47\text{ Å}$). Mettendo il protone su `ND1` (**`HSD`**), `HIS 40` dona un ponte a idrogeno inter-repeat fondamentale a `SER 74` che tiene cuciti insieme `ANK2` e `ANK3`!
     - **`HIS 67` (nel $\beta$-hairpin loop tra `ANK2` e `ANK3`) $\rightarrow$ `HSE`:** nel file NMR originale X-PLOR aveva messo di default il protone `HD1` su `ND1`. Ma a **`3.24 Å` da `ND1`** si trova il gruppo ammonico carico positivamente **`NZ` ($-\text{NH}_3^+$) di `LYS 37`** (che può solo *donare* H-bond!), mentre a **`2.40 Å` da `NE2`** si trova l'ossigeno carbonilico **`O` (`C=O`) di `GLY 35`** (che può solo *accettare* H-bond!). Spostando il protone su `NE2` (**`HSE`**), abbiamo eliminato una repulsione `+`/`+` artificiale tra `LYS 37` e `HIS 67` e attivato una **doppia cucitura cooperativa a idrogeno**: `LYS 37 (NZ) -> HIS 67 (ND1)` e `HIS 67 (NE2-HE2) -> GLY 35 (O)`!
     - **`HIS 68` e `HIS 80` $\rightarrow$ `HSD`:** esposte al solvente o in contatto apolare (`LEU 48`), mantengono lo stato `HSD` nativo dell'NMR.

### 4.4 La "Functional Mask" Oggettiva (24 Residui Intoccabili)
Dal cutoff atomico $< 4.5\text{ Å}$ nel complesso `3AAA`, abbiamo estratto i **24 residui funzionali** che costituiscono l'impronta di legame sulla faccia concava della Myotrophin:
- **`ANK1` (12 residui):** `LYS 4, GLU 5, MET 7, TRP 8, ALA 9, LYS 11, ASN 12, ASP 14, ASP 16, GLU 17, GLY 34, ARG 36`
- **`ANK2` (7 residui):** `TYR 41, ASP 44, CYS 45, GLN 47, LYS 66, HIS 67, ILE 69`
- **`ANK3` (5 residui):** `TYR 77, GLU 78, HIS 80, PRO 99, ASP 100`
- **`ANK4` (0 residui):** Il C-cap (`101–118`) non tocca l'Actin-Capping Protein: ha una funzione puramente strutturale di nucleo di folding!
- **Considerazione personale:** Questa maschera include al 100% la triade calda sperimentale (**`TRP 8, ASP 44, GLU 78`**) dimostrata da *Takeda et al. (2010)* tramite mutagenesi sito-specifica, e ci lascia esattamente **94 residui non-funzionali** (88 interni + 6 terminali) su cui ingegnerizzare la termostabilità con la certezza matematica di non alterare l'affinità per l'Actin-Capping Protein.

---

## 5. ⏭️ Come Agire Dopo nel Futuro (Impatto sui Punti Successivi e sulla Presentazione)

1. **Impatto diretto sul Punto 2 (Setup MD):** I tautomeri verificati (`HIS 40 -> HSD`, `HIS 67 -> HSE`, `HIS 68/80 -> HSD`) e la rimozione degli idrogeni NMR (`protein and noh`) vengono codificati direttamente in [`scripts/02_build_wt_charmm36m.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/02_build_wt_charmm36m.tcl) (e verranno riapplicati in modo identico a tutti i mutanti nel Punto 4 per garantire la *Matched Comparison*).
2. **Impatto diretto sul Punto 3 (Selezione Mutazioni):**
   - La lista dei 24 residui della Functional Mask viene caricata come filtro automatico (*Regola 1*) in tutti gli script di screening ([`04_ankyrin_consensus_and_targets.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/04_ankyrin_consensus_and_targets.py), [`05_parse_foldx_screening.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/05_parse_foldx_screening.py), [`06_build_consensus_matrix.py`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/06_build_consensus_matrix.py)), permettendoci di scartare subito false proposte di FireProt come `D16W` o `I69G/M/R`.
   - Inoltre, la scoperta del ponte a idrogeno inter-repeat **`HIS 40 (HSD) -> SER 74 (OG)`** ci ha permesso nel Punto 3 di riconoscere e scartare la "trappola" **`S74A`** (che avrebbe tolto l'ossigeno accettore `OG` a `SER 74`)!
3. **Cosa portare nelle Slide della Presentazione del 16 Dicembre (Slide 1 & 2, ~2 minuti):**
   - Rendering 3D VMD da [`scripts/viz_01_structure_and_mask.tcl`](file:///home/olly/UNI/BioMolEng2026/laboratorio/BML/scripts/viz_01_structure_and_mask.tcl) che mostra i 4 Ankyrin Repeats colorati (`ANK1` Blu, `ANK2` Ciano, `ANK3` Arancio, `ANK4` Verde), la nuvola rossa semitrasparente della **Functional Mask** (`3AAA` $< 4.5\text{ Å}$ con la triade `Trp8, Asp44, Glu78`) e lo zoom chimico sul doppio ponte a idrogeno di `HIS 40 (HSD)` e `HIS 67 (HSE)`.
