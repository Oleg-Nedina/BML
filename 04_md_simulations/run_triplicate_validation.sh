#!/usr/bin/env bash
# ==============================================================================
# run_triplicate_validation.sh
# Script batch per il TRIPLICATO MULTI-SEED (n = 3 velocity seeds a 400 K, 2.5 ns):
#   - Sfrutta il Seed 1 già completato come Replica 1 su 3 (protocollo identico
#     confermato dal Prof. Gautieri il 25 Nov 2026!).
#   - Esegue i restanti Seed 2 (23456) e Seed 3 (34567) con resume automatico.
#
# Modalità disponibili:
#   1) wt_triple   -> 4 run nuove: WT (s2, s3) + Triplo K37L_H68G_T103L (s2, s3)  [CONSIGLIATA!]
#   2) core_triple -> 6 run nuove: WT (s2, s3) + Doppio (s2, s3) + Triplo (s2, s3)
#   3) oleg        -> 2 run nuove: Triplo K37L_H68G_T103L (s2, s3)   [1 notte su PC Oleg]
#   4) alessia     -> 2 run nuove: WT (s2, s3)                        [1 notte su PC Alessia]
#   5) all         -> 10 run nuove: s2 e s3 per tutti e 5 i sistemi
#
# Al termine lancia automaticamente:
#   - scripts/09_analyze_matched_validation.tcl  (Q-value Yang 1998, RMSD, RMSF)
#   - scripts/12_analyze_hbonds_and_helices.tcl  (H-Bonds 3.5A/30° + 7 Alfa-Eliche H1-H7)
#   - scripts/10_plot_final_validation_bekker.py + 11 + 12 (Tabelle μ ± σ + Grafici PDF/PNG)
#
# Uso dalla cartella BML/:
#   bash 04_md_simulations/run_triplicate_validation.sh wt_triple 8
# ==============================================================================

set -euo pipefail

MODE="${1:-wt_triple}"
NCORES="${2:-8}"
NAMD3_BIN="${NAMD3_BIN:-$(command -v namd3 || echo /usr/local/bin/namd3)}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BML_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$BML_ROOT"

# Normalizza eventuali percorsi relativi nei file .conf
for cf in 04_md_simulations/*/*.conf; do
    [[ -f "$cf" ]] || continue
    sed -i -E \
        -e 's|([[:space:]])04_md_simulations/[^/]+/|\1|g' \
        -e 's|([[:space:]])01_structures/|\1../../01_structures/|g' \
        -e 's|([[:space:]])common/charmm36/|\1../../common/charmm36/|g' \
        "$cf"
done

case "$MODE" in
    wt_triple)
        CONFS=(
            "04_md_simulations/WT/wt_400K_seed2.conf"
            "04_md_simulations/WT/wt_400K_seed3.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed3.conf"
        )
        ;;
    core_triple)
        CONFS=(
            "04_md_simulations/WT/wt_400K_seed2.conf"
            "04_md_simulations/WT/wt_400K_seed3.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed3.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed3.conf"
        )
        ;;
    oleg)
        CONFS=(
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed3.conf"
        )
        ;;
    alessia)
        CONFS=(
            "04_md_simulations/WT/wt_400K_seed2.conf"
            "04_md_simulations/WT/wt_400K_seed3.conf"
        )
        ;;
    all)
        CONFS=(
            "04_md_simulations/WT/wt_400K_seed2.conf"
            "04_md_simulations/WT/wt_400K_seed3.conf"
            "04_md_simulations/K37L/K37L_400K_seed2.conf"
            "04_md_simulations/K37L/K37L_400K_seed3.conf"
            "04_md_simulations/T103L/T103L_400K_seed2.conf"
            "04_md_simulations/T103L/T103L_400K_seed3.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed3.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed3.conf"
        )
        ;;
    *)
        echo "Modalita' non valida: $MODE (usa: wt_triple | core_triple | oleg | alessia | all)"
        exit 1
        ;;
esac

echo "=========================================================================="
echo " [$(date)] AVVIO TRIPLICATO MULTI-SEED — Modalita': $MODE (namd3 +p${NCORES})"
echo " Configurazioni in coda: ${#CONFS[@]}"
echo "=========================================================================="

idx=0
total="${#CONFS[@]}"
for conf in "${CONFS[@]}"; do
    idx=$((idx + 1))
    base="${conf%.conf}"
    dcd="${base}.dcd"
    log="${base}.log"

    if [[ -f "$dcd" && -f "$log" ]] && grep -q "End of program" "$log" 2>/dev/null && ! grep -q "FATAL ERROR" "$log" 2>/dev/null; then
        echo " -> [$idx/$total] [SKIP] Già completata: $dcd"
        continue
    fi

    echo " -> [$idx/$total] [$(date)] Avvio NAMD 3: $conf ..."
    "$NAMD3_BIN" +p"${NCORES}" "$conf" > "$log" 2>&1
    echo "    [$(date)] Completata: $dcd !"
done

echo "=========================================================================="
echo " [$(date)] Aggiornamento Analisi VMD (Q-value Yang 1998 + H-Bonds + 7 Eliche)..."
echo "=========================================================================="
vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl > 05_analysis/vmd_triplicate_q.log 2>&1
vmd -dispdev text -e scripts/12_analyze_hbonds_and_helices.tcl > 05_analysis/vmd_triplicate_helices.log 2>&1

python3 scripts/10_plot_final_validation_bekker.py
python3 scripts/11_compare_wt_vs_triple_for_prof.py
python3 scripts/12_plot_hbonds_and_helices.py

echo "=========================================================================="
echo " [$(date)] TRIPLICATO '$MODE' E TUTTI I GRAFICI COMPLETATI!"
echo "=========================================================================="
