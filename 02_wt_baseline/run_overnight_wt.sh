#!/usr/bin/env bash
# ==============================================================================
# run_overnight_wt.sh
# Script batch per eseguire durante la notte l'intera Fase 2 (ROAD2.md):
#   1. Simulazione MD Wild-Type a 300 K (1.5 ns, controllo fisiologico)
#   2. Simulazione MD Wild-Type a 400 K (2.5 ns, stress termico Bekker et al.)
#   3. Analisi automatica in VMD (RMSD, RMSF, Q-value globale/chimico/per-repeat, Q_i)
#   4. Generazione automatica dei grafici in Python per le slide (06_presentation/)
#
# Esecuzione dalla cartella principale BML/:
#   bash 02_wt_baseline/run_overnight_wt.sh
# (Oppure in background resistente alla chiusura del terminale:
#   nohup bash 02_wt_baseline/run_overnight_wt.sh > 02_wt_baseline/overnight_master.log 2>&1 & )
# ==============================================================================

set -e

# Imposta il numero di thread (default 8 core fisici, +40% piu' veloce di 16 SMT su AVX-512!)
NPROCS=${1:-8}

echo "=================================================================="
echo " [$(date)] AVVIO BATCH NOTTURNO FASE 2 (NAMD3 +p${NPROCS})"
echo "=================================================================="

echo " -> [1/4] Avvio WT a 300 K (750.000 step = 1.5 ns)..."
namd3 +p${NPROCS} 02_wt_baseline/wt_300K.conf > 02_wt_baseline/wt_300K.log
echo "    [$(date)] Completata run WT 300 K!"

echo " -> [2/4] Avvio WT a 400 K (1.250.000 step = 2.5 ns)..."
namd3 +p${NPROCS} 02_wt_baseline/wt_400K.conf > 02_wt_baseline/wt_400K.log
echo "    [$(date)] Completata run WT 400 K!"

echo " -> [3/4] Esecuzione analisi traiettorie in VMD (03_analyze_wt_baseline.tcl)..."
vmd -dispdev text -e scripts/03_analyze_wt_baseline.tcl > 02_wt_baseline/analysis_wt.log
echo "    [$(date)] Analisi VMD completata! Report salvato in 02_wt_baseline/unstable_residues_report.txt"

echo " -> [4/4] Generazione grafici Python per le slide (plot_wt_baseline.py)..."
python3 scripts/plot_wt_baseline.py
echo "=================================================================="
echo " [$(date)] TUTTO COMPLETATO CON SUCCESSO!"
echo " Controlla:"
echo "   - 02_wt_baseline/unstable_residues_report.txt          (Classifica bulloni deboli)"
echo "   - 06_presentation/fig5_wt_baseline_plots.pdf           (Grafici vettoriali Slide 4)"
echo "   - 06_presentation/fig5b_wt_analysisMD_observables.pdf  (Grafici analysisMD.tcl)"
echo "=================================================================="
