#!/usr/bin/env bash
# ==============================================================================
# run_4_screening_400K.sh
# Esegue in sequenza le 4 simulazioni di screening a 400 K (Seed 1, 2.5 ns):
#   1. Simulazione 1: K37L             (Mutante Singolo 1, ANK2)
#   2. Simulazione 2: T103L            (Mutante Singolo 2, ANK4)
#   3. Simulazione 3: K37L_T103L       (Doppio Mutante, ANK2 + ANK4)
#   4. Simulazione 4: K37L_H68G_T103L  (Triplo Mutante, ANK2 + ANK3 + ANK4)
#
# Funzionalità integrate:
#   - Prima di iniziare, analizza e genera già i 3 plot + report separati del
#     Wild-Type (05_analysis/single_runs/00_WT/).
#   - Salva tutti i file di output NAMD (.dcd, .log, .xst, .restart.*)
#     ordinatamente nella cartella di ciascun mutante (04_md_simulations/<MUT>/).
#   - Resume automatico: se una delle 4 simulazioni è già stata completata
#     ("End of program" nel .log e file .dcd presente), la salta senza rifarla.
#   - Subito dopo il completamento di OGNI singola simulazione, lancia VMD
#     (09_analyze_matched_validation.tcl) e Python (10_plot_single_runs_400K.py)
#     per generare e salvare immediatamente i 3 plot PDF individuali e il report
#     di quella simulazione in 05_analysis/single_runs/<NUMERO_SIMULAZIONE>/!
#
# Uso dalla cartella BML/:
#   ./04_md_simulations/run_4_screening_400K.sh
# ==============================================================================

set -euo pipefail

NAMD3_BIN="${NAMD3_BIN:-$(command -v namd3 || echo /usr/local/bin/namd3)}"
NAMD_FLAGS="+p8"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BML_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$BML_ROOT"

STATUS_LOG="$BML_ROOT/04_md_simulations/screening_4_runs_progress.log"

log_msg() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "$msg" | tee -a "$STATUS_LOG"
}

update_single_plots() {
    log_msg ">>> Aggiornamento analisi VMD e generazione plot singoli separati..."
    vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl > 05_analysis/vmd_single_analysis.log 2>&1
    python3 scripts/10_plot_single_runs_400K.py | tee -a "$STATUS_LOG"
}

log_msg "=========================================================================="
log_msg " AVVIO PIPELINE 4 SIMULAZIONI SCREENING A 400 K (SEED 1, 2.5 ns CIASCUNA)"
log_msg "=========================================================================="

# Normalizza tutti i percorsi relativi nei file .conf (poiché NAMD fa chdir in 04_md_simulations/<VAR>/)
for cf in "$BML_ROOT"/04_md_simulations/*/*.conf; do
    [[ -f "$cf" ]] || continue
    sed -i -E \
        -e 's|([[:space:]])04_md_simulations/[^/]+/|\1|g' \
        -e 's|([[:space:]])01_structures/|\1../../01_structures/|g' \
        -e 's|([[:space:]])common/charmm36/|\1../../common/charmm36/|g' \
        "$cf"
done

# Passo 0: Genera subito i plot individuali del Wild-Type (02_wt_baseline/wt_400K.dcd)
log_msg ">>> Passo 0: Verifica e generazione plot individuali per il Wild-Type (00_WT)..."
update_single_plots

# Lista ordinata delle 4 simulazioni da eseguire
SIM_NAMES=(
    "K37L"
    "T103L"
    "K37L_T103L"
    "K37L_H68G_T103L"
)

SIM_LABELS=(
    "Simulazione 1/4: Mutante Singolo K37L (ANK2)"
    "Simulazione 2/4: Mutante Singolo T103L (ANK4)"
    "Simulazione 3/4: Doppio Mutante K37L_T103L (ANK2 + ANK4)"
    "Simulazione 4/4: Triplo Mutante K37L_H68G_T103L (ANK2 + ANK3 + ANK4)"
)

for idx in "${!SIM_NAMES[@]}"; do
    sys="${SIM_NAMES[$idx]}"
    label="${SIM_LABELS[$idx]}"
    wdir="$BML_ROOT/04_md_simulations/$sys"
    conf_f="04_md_simulations/${sys}/${sys}_400K_seed1.conf"
    log_f="04_md_simulations/${sys}/${sys}_400K_seed1.log"
    dcd_f="04_md_simulations/${sys}/${sys}_400K_seed1.dcd"

    log_msg "--------------------------------------------------------------------------"
    log_msg " $label"
    log_msg " Config   : $conf_f"
    log_msg "--------------------------------------------------------------------------"

    if [[ -f "$BML_ROOT/$dcd_f" && -f "$BML_ROOT/$log_f" ]] && grep -q "End of program" "$BML_ROOT/$log_f" 2>/dev/null && ! grep -q "FATAL ERROR" "$BML_ROOT/$log_f" 2>/dev/null; then
        log_msg "[SKIP] $sys (Seed 1, 400K) è già completata ($dcd_f presente con 'End of program')."
    else
        rm -f "$wdir/${sys}_400K_s1.log"
        log_msg "[START] Esecuzione NAMD 3 per $sys ($conf_f)..."
        "$NAMD3_BIN" $NAMD_FLAGS "$conf_f" > "$log_f" 2>&1
        if [[ -f "$BML_ROOT/$dcd_f" ]] && grep -q "End of program" "$BML_ROOT/$log_f" 2>/dev/null && ! grep -q "FATAL ERROR" "$BML_ROOT/$log_f" 2>/dev/null; then
            speed_info=$(grep "Info: Benchmark time:" "$BML_ROOT/$log_f" | tail -n 1 || true)
            log_msg "[DONE] $sys completata con successo! $speed_info"
        else
            log_msg "[ERRORE] $sys non ha terminato correttamente. Controlla $log_f"
            exit 1
        fi
    fi

    # Genera subito i plot singoli e il report per la simulazione appena conclusa
    update_single_plots
done

log_msg "=========================================================================="
log_msg " TUTTE E 4 LE SIMULAZIONI E I RELATIVI PLOT SINGOLI SONO STATI COMPLETATI!"
log_msg " Cartelle con i plot separati:"
log_msg "   - 05_analysis/single_runs/00_WT/"
log_msg "   - 05_analysis/single_runs/01_K37L/"
log_msg "   - 05_analysis/single_runs/02_T103L/"
log_msg "   - 05_analysis/single_runs/03_K37L_T103L/"
log_msg "   - 05_analysis/single_runs/04_K37L_H68G_T103L/"
log_msg "=========================================================================="
