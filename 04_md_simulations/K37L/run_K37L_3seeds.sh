#!/usr/bin/env bash
# ==============================================================================
# run_K37L_3seeds.sh  (Sistema: Singola Mutazione 1 — K37L in Repeat 1)
# Esegue i 3 velocity seeds identici a 400 K (2.5 ns ciascuno) per l'analisi di varianza:
#   - Seed 1 (12345): K37L_400K_seed1.conf  [già completato -> SKIP automatico]
#   - Seed 2 (23456): K37L_400K_seed2.conf
#   - Seed 3 (34567): K37L_400K_seed3.conf
#
# Uso dalla cartella BML/:
#   bash 04_md_simulations/K37L/run_K37L_3seeds.sh 8
# ==============================================================================
set -euo pipefail

NCORES="${1:-8}"
NAMD3_BIN="${NAMD3_BIN:-$(command -v namd3 || echo /usr/local/bin/namd3)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BML_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$BML_ROOT"

SYS_DIR="04_md_simulations/K37L"
PREFIX="K37L_400K"

for cf in "$SYS_DIR"/${PREFIX}_seed*.conf; do
    [[ -f "$cf" ]] || continue
    sed -i -E \
        -e 's|([[:space:]])04_md_simulations/[^/]+/|\1|g' \
        -e 's|([[:space:]])01_structures/|\1../../01_structures/|g' \
        -e 's|([[:space:]])common/charmm36/|\1../../common/charmm36/|g' \
        "$cf"
done

echo "=========================================================================="
echo " [$(date)] TRIPLICATO 400 K — K37L (3 seeds: 12345, 23456, 34567)"
echo "=========================================================================="

for s in 1 2 3; do
    conf="${SYS_DIR}/${PREFIX}_seed${s}.conf"
    dcd="${SYS_DIR}/${PREFIX}_seed${s}.dcd"
    log="${SYS_DIR}/${PREFIX}_seed${s}.log"

    if [[ -f "$dcd" && -f "$log" ]] && grep -q "End of program" "$log" 2>/dev/null && ! grep -q "FATAL ERROR" "$log" 2>/dev/null; then
        echo " -> [Seed $s/3] [SKIP] Già completata: $dcd"
        continue
    fi

    echo " -> [Seed $s/3] [$(date)] Avvio NAMD 3 (+p${NCORES}): $conf ..."
    "$NAMD3_BIN" +p"${NCORES}" "$conf" > "$log" 2>&1
    echo "    [$(date)] Completata: $dcd !"
done

echo "=========================================================================="
echo " [$(date)] TUTTI E 3 I SEED DI K37L SONO COMPLETATI!"
echo "=========================================================================="
