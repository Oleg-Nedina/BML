#!/usr/bin/env bash
# ==============================================================================
# 04b_run_remaining_foldx_parallel.sh
# Scansiona con FoldX PositionScan le restanti 45 posizioni interne dello scaffold
# (portando la copertura da 43/88 al 100% = 88/88 posizioni modificabili!),
# dividendo il lavoro su 8 core in parallelo (~3-4 minuti totali invece di 28 min!).
#
# Esecuzione dalla cartella principale BML/:
#   bash scripts/04b_run_remaining_foldx_parallel.sh
# ==============================================================================

set -e

echo "=================================================================="
echo " AVVIO SCANSIONE PARALLELA (8 CORE) DELLE 45 POSIZIONI RIMANENTI"
echo "=================================================================="

# 1. Pulisci prima i ~1000 PDB temporanei della prima run
rm -f 03_mutations/[A-Z]*[0-9]*_2MYO_Repair.pdb

# 2. Genera le 8 sottocartelle di lavoro isolate per non far confliggere FoldX
python3 - << 'EOF'
import os

SEQ_2MYO = (
    "MCDKEFMWALKNGDLDEVKDYVAKGEDVNRTLE"
    "GGRKPLHYAADCGQLEILEFLLLKGADINAPDK"
    "HHITPLLSAVYEGHVSCVKLLLSKGADKTVKGP"
    "DGLTALEATDNQAIKALLQ"
)
FUNCTIONAL_MASK = {
    4, 5, 7, 8, 9, 11, 12, 14, 16, 17, 34, 36,
    41, 44, 45, 47, 66, 67, 69, 77, 78, 80, 99, 100
}
ALREADY_SCANNED = {
    24, 26, 27, 28, 29, 31, 32, 33, 35, 37, 38, 40, 48, 53, 54,
    57, 59, 60, 62, 64, 68, 70, 71, 72, 73, 81, 84, 85, 87, 88,
    91, 92, 93, 94, 95, 96, 97, 101, 103, 104, 105, 108, 110
}

remaining = sorted(list(set(range(4, 116)) - FUNCTIONAL_MASK - ALREADY_SCANNED))
print(f"Totale posizioni interne rimanenti da scansionare: {len(remaining)} -> {remaining}")

n_workers = 8
chunks = [remaining[i::n_workers] for i in range(n_workers)]

for idx, chunk in enumerate(chunks, 1):
    wdir = f"03_mutations/worker_{idx}"
    os.makedirs(wdir, exist_ok=True)
    # Link o copia 2MYO_Repair.pdb
    target_pdb = os.path.join(wdir, "2MYO_Repair.pdb")
    if not os.path.exists(target_pdb):
        os.system(f"cp 03_mutations/2MYO_Repair.pdb {target_pdb}")
    codes = [f"{SEQ_2MYO[r-1]}A{r}d" for r in chunk]
    with open(os.path.join(wdir, "scan.cfg"), "w") as f:
        f.write("command=PositionScan\n")
        f.write("pdb-dir=.\n")
        f.write("output-dir=.\n")
        f.write("pdb=2MYO_Repair.pdb\n")
        f.write("positions=" + ",".join(codes) + "\n")
    print(f"  Worker {idx}: {len(chunk)} posizioni ({','.join(codes)})")
EOF

# 3. Lancia gli 8 processi FoldX isolati dentro le rispettive worker_* (~3-10 minuti)
echo " -> Esecuzione in parallelo su 8 core..."
pids=()
for i in {1..8}; do
    (cd "03_mutations/worker_${i}" && foldx -f scan.cfg > scan.log 2>&1) &
    pids+=($!)
done

# Attendi il completamento di tutti gli 8 worker
for pid in "${pids[@]}"; do
    wait "$pid"
done
echo " -> Tutti gli 8 worker FoldX hanno terminato!"

# 4. Unisci i risultati nel file principale e copia i file energies_*
for i in {1..8}; do
    if [ -f "03_mutations/worker_${i}/PS_2MYO_Repair_scanning_output.txt" ]; then
        cat "03_mutations/worker_${i}/PS_2MYO_Repair_scanning_output.txt" >> 03_mutations/PS_2MYO_Repair_scanning_output.txt
    fi
    cp -f 03_mutations/worker_${i}/energies_*_2MYO_Repair.txt 03_mutations/ 2>/dev/null || true
    rm -rf "03_mutations/worker_${i}"
done

# 5. Aggiorna la classifica globale (88/88 posizioni!) e la Matrice di Consenso
python3 scripts/05_parse_foldx_screening.py
python3 scripts/06_build_consensus_matrix.py

echo "=================================================================="
echo " COPERTURA 100% COMPLETATA (88/88 posizioni interne valutate)!"
echo " Controlla: 03_mutations/consensus_master_table.txt"
echo "=================================================================="
