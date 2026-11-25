#!/usr/bin/env bash
# ==============================================================================
# genera_confronto_prof_WT_vs_Triplo.sh
# Crea/aggiorna in un solo comando la cartella pronta per il colloquio col Prof:
#   06_presentation/presentazione_prof_WT_vs_Triplo/
#
# Esegue:
#   1. Verifica/aggiorna l'analisi VMD con riferimento Frame 0 = DCD (Q(0) = 1.00)
#   2. Genera la tabella 1-vs-1 e i 3 grafici comparativi WT vs Triplo (.pdf + .png)
#   3. Copia nella stessa cartella anche il grafico FoldX BuildModel e la traccia
#      spiegazione_sintetica.md così hai tutto in una sola directory!
#
# Uso dalla cartella BML/:
#   bash scripts/genera_confronto_prof_WT_vs_Triplo.sh
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BML_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$BML_ROOT"

OUT_DIR="06_presentation/presentazione_prof_WT_vs_Triplo"
mkdir -p "$OUT_DIR"

echo "=========================================================================="
echo " 1/3 — Verifica analisi VMD (Frame 0 = 1° frame DCD, Q(0) = 1.00)..."
echo "=========================================================================="
vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl

echo "=========================================================================="
echo " 2/3 — Generazione plot e tabella comparativa WT vs Triplo Mutante..."
echo "=========================================================================="
python3 scripts/11_compare_wt_vs_triple_for_prof.py

echo "=========================================================================="
echo " 3/3 — Copia file di supporto (FoldX BuildModel + Discorso) nella cartella"
echo "=========================================================================="
[[ -f "06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf" ]] && \
    cp "06_presentation/fig8b_buildmodel_singles_doubles_epistasis.pdf" "$OUT_DIR/04_FoldX_BuildModel_singoli_doppi_triplo.pdf"
[[ -f "spiegazione_sintetica.md" ]] && \
    cp "spiegazione_sintetica.md" "$OUT_DIR/spiegazione_sintetica.md"

echo ""
cat "$OUT_DIR/00_TABELLA_CONFRONTO_WT_vs_TRIPLO.txt"
echo ""
echo ">>> TUTTO PRONTO! Apri la cartella: $OUT_DIR/"
ls -lh "$OUT_DIR/"
