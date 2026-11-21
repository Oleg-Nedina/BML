#!/usr/bin/env bash
# ==============================================================================
# 04_md_simulations/run_batch_phase4.sh
# Script Batch Intelligente per le Simulazioni MD Multi-Seed della Fase 4 (400 K)
#
# Caratteristiche:
#   - Salta automaticamente qualsiasi simulazione il cui file .dcd esista già
#     e sia completo (così non ricalcola mai una run già finita!).
#   - Supporta 5 modalità di esecuzione per dividere comodamente il lavoro:
#       1) "seed1"   -> Esegue il Seed 1 (12345) dei 3 mutanti principali:
#                       K37L, T103L, K37L_T103L (3 run -> confronto immediato col WT Seed 1!)
#       2) "oleg"    -> Esegue i Seed 2 e 3 di K37L e K37L_T103L (4 run per il PC Linux di Oleg)
#       3) "alessia" -> Esegue i Seed 2 e 3 di WT e T103L (4 run per il PC Windows di Alessia)
#       4) "all"     -> Esegue tutte le 11 run del quartetto principale (WT, K37L, T103L, K37L_T103L)
#       5) "bonus"   -> Esegue il Seed 1 di H68G (#1 BTC FireProt) e del Triplo K37L_H68G_T103L
#   - Al termine delle simulazioni, lancia automaticamente l'analisi VMD
#     (scripts/09_analyze_matched_validation.tcl) e genera la Table II di Bekker
#     e i grafici PDF (scripts/10_plot_final_validation_bekker.py).
#
# Uso dalla cartella principale BML/:
#   nohup bash 04_md_simulations/run_batch_phase4.sh seed1 8 > 04_md_simulations/phase4_seed1.log 2>&1 &
# ==============================================================================

set -e

MODE="${1:-seed1}"
NCORES="${2:-8}"

case "$MODE" in
    seed1)
        CONFS=(
            "04_md_simulations/K37L/K37L_400K_seed1.conf"
            "04_md_simulations/T103L/T103L_400K_seed1.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed1.conf"
        )
        ;;
    oleg)
        CONFS=(
            "04_md_simulations/K37L/K37L_400K_seed1.conf"
            "04_md_simulations/T103L/T103L_400K_seed1.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed1.conf"
            "04_md_simulations/K37L/K37L_400K_seed2.conf"
            "04_md_simulations/K37L/K37L_400K_seed3.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed3.conf"
        )
        ;;
    alessia)
        CONFS=(
            "04_md_simulations/WT/wt_400K_seed2.conf"
            "04_md_simulations/WT/wt_400K_seed3.conf"
            "04_md_simulations/T103L/T103L_400K_seed2.conf"
            "04_md_simulations/T103L/T103L_400K_seed3.conf"
        )
        ;;
    all)
        CONFS=(
            "04_md_simulations/K37L/K37L_400K_seed1.conf"
            "04_md_simulations/T103L/T103L_400K_seed1.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed1.conf"
            "04_md_simulations/WT/wt_400K_seed2.conf"
            "04_md_simulations/WT/wt_400K_seed3.conf"
            "04_md_simulations/K37L/K37L_400K_seed2.conf"
            "04_md_simulations/K37L/K37L_400K_seed3.conf"
            "04_md_simulations/T103L/T103L_400K_seed2.conf"
            "04_md_simulations/T103L/T103L_400K_seed3.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed2.conf"
            "04_md_simulations/K37L_T103L/K37L_T103L_400K_seed3.conf"
        )
        ;;
    bonus)
        CONFS=(
            "04_md_simulations/H68G/H68G_400K_seed1.conf"
            "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed1.conf"
        )
        ;;
    *)
        echo "Modalita' non riconosciuta: $MODE (usa: seed1 | oleg | alessia | all | bonus)"
        exit 1
        ;;
esac

echo "=========================================================================="
echo " [$(date)] AVVIO BATCH FASE 4 — Modalita': $MODE (namd3 +p${NCORES})"
echo " Totale configurazioni nella coda: ${#CONFS[@]}"
echo "=========================================================================="

idx=0
total="${#CONFS[@]}"
for conf in "${CONFS[@]}"; do
    idx=$((idx + 1))
    base="${conf%.conf}"
    dcd="${base}.dcd"
    log="${base}.log"

    if [ -f "$dcd" ] && [ -f "$log" ] && grep -q "End of program" "$log" 2>/dev/null; then
        echo " -> [$idx/$total] GIÀ COMPLETATA: $dcd (salto)"
        continue
    fi

    echo " -> [$idx/$total] [$(date)] Avvio NAMD3: $conf ..."
    namd3 +p"${NCORES}" "$conf" > "$log" 2>&1
    echo "    [$(date)] Completata: $dcd !"
done

echo "=========================================================================="
echo " [$(date)] Esecuzione Analisi VMD Multi-Seed (09_analyze_matched_validation.tcl)..."
echo "=========================================================================="
vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl > 05_analysis/analysis_matched.log 2>&1

echo "=========================================================================="
echo " [$(date)] Generazione Bekker Table II e Grafici PDF (10_plot_final_validation_bekker.py)..."
echo "=========================================================================="
python3 scripts/10_plot_final_validation_bekker.py

echo "=========================================================================="
echo " [$(date)] BATCH '$MODE' TERMINATO CON SUCCESSO!"
echo " Controlla: cat 05_analysis/bekker_table2_and_epistasis_report.txt"
echo "=========================================================================="
