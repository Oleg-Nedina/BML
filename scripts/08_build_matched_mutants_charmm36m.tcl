# ==============================================================================
# 08_build_matched_mutants_charmm36m.tcl
# Punto 4 della Scaletta del Prof. Gautieri (Laboratory_Project_2026.pdf, Slide 3, 15, 21)
# & Livello 2 della Politica di Validazione (ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md):
#
# Costruisce per tutte le varianti selezionate (Mut1 = K37L, Mut2 = T103L,
# Double = K37L_T103L, piu' opzionali H68G e Triple = K37L_H68G_T103L)
# esattamente lo STESSO dominio solvatato e ionizzato del Wild-Type ("Matched Protocol"):
#   1. Assegnazione tautomerica identica delle Istidine verificata in Fase 1:
#      - HIS 40 -> HSD (ponte H inter-repeat con SER 74 OG)
#      - HIS 67 -> HSE (doppio ponte H con LYS/LEU 37 e GLY 35 O)
#      - HIS 68 -> HSD (se non mutata in GLY 68)
#      - HIS 80 -> HSD
#   2. Rimozione idrogeni di FoldX ("protein and noh") e ricostruzione con psfgen
#      usando la topologia CHARMM36m (common/charmm36/top_all36_prot.rtf).
#   3. Solvatazione esplicita TIP3P con padding 12.0 A (solvate -t 12).
#   4. Neutralizzazione elettrica e forza ionica 0.15 M NaCl (autoionize -sc 0.15).
#   5. Generazione automatica dei file di configurazione NAMD (.conf) per i 3 seed
#      indipendenti a 400 K (seed 12345, 23456, 34567) in 04_md_simulations/<variant>/,
#      piu' i 2 seed rimanenti del Wild-Type (seed 23456, 34567).
#
# Esecuzione dalla cartella principale BML/:
#   vmd -dispdev text -e scripts/08_build_matched_mutants_charmm36m.tcl
# ==============================================================================

package require psfgen
package require solvate
package require autoionize

file mkdir "04_md_simulations"

# Procedura helper per generare un file .conf di NAMD 3 a 400 K perfettamente "Matched" con wt_400K.conf
proc write_matched_namd_conf {conf_path psf_rel pdb_rel out_rel seed_val lx ly lz cx cy cz} {
    set fp [open $conf_path w]
    puts $fp "# =============================================================================="
    puts $fp "# Configurazione NAMD 3 Matched Comparison (400 K, 2.5 ns = 1.250.000 step)"
    puts $fp "# Generata automaticamente da scripts/08_build_matched_mutants_charmm36m.tcl"
    puts $fp "# Protocollo identico a 02_wt_baseline/wt_400K.conf (Slide 3, 15, 21)"
    puts $fp "# =============================================================================="
    puts $fp "structure          $psf_rel"
    puts $fp "coordinates        $pdb_rel"
    puts $fp "set temperature    400"
    puts $fp "set outputname     $out_rel"
    puts $fp "firsttimestep      0"
    puts $fp ""
    puts $fp "paraTypeCharmm     on"
    puts $fp "parameters         common/charmm36/par_all36m_prot.prm"
    puts $fp "parameters         common/charmm36/toppar_water_ions_namd.prm"
    puts $fp "temperature        \$temperature"
    puts $fp "seed               $seed_val"
    puts $fp ""
    puts $fp "# Force-Field Parameters (Matched)"
    puts $fp "exclude            scaled1-4"
    puts $fp "1-4scaling         1.0"
    puts $fp "cutoff             12.0"
    puts $fp "switching          on"
    puts $fp "switchdist         10.0"
    puts $fp "pairlistdist       14.0"
    puts $fp ""
    puts $fp "# Integrator Parameters (dt = 2.0 fs, SHAKE/RATTLE)"
    puts $fp "timestep           2.0"
    puts $fp "rigidBonds         all"
    puts $fp "nonbondedFreq      1"
    puts $fp "fullElectFrequency 2"
    puts $fp "stepspercycle      10"
    puts $fp ""
    puts $fp "# Constant Temperature Control (Langevin Thermostat)"
    puts $fp "langevin           on"
    puts $fp "langevinDamping    1"
    puts $fp "langevinTemp       \$temperature"
    puts $fp "langevinHydrogen   off"
    puts $fp ""
    puts $fp "# Periodic Boundary Conditions"
    puts $fp [format "cellBasisVector1   %.3f   0.000   0.000" $lx]
    puts $fp [format "cellBasisVector2   0.000   %.3f   0.000" $ly]
    puts $fp [format "cellBasisVector3   0.000   0.000   %.3f" $lz]
    puts $fp [format "cellOrigin         %.3f   %.3f   %.3f" $cx $cy $cz]
    puts $fp "wrapAll            on"
    puts $fp ""
    puts $fp "# PME (Full-System Periodic Electrostatics)"
    puts $fp "PME                yes"
    puts $fp "PMEGridSpacing     1.0"
    puts $fp ""
    puts $fp "# Constant Pressure Control (Nose-Hoover Langevin Piston, 1 atm)"
    puts $fp "useGroupPressure      yes"
    puts $fp "useFlexibleCell       no"
    puts $fp "useConstantArea       no"
    puts $fp "langevinPiston        on"
    puts $fp "langevinPistonTarget  1.01325"
    puts $fp "langevinPistonPeriod  100.0"
    puts $fp "langevinPistonDecay   50.0"
    puts $fp "langevinPistonTemp    \$temperature"
    puts $fp ""
    puts $fp "# Output Frequencies (1 frame ogni 1000 step = 2 ps -> 1250 frame in 2.5 ns)"
    puts $fp "outputName         \$outputname"
    puts $fp "restartfreq        5000"
    puts $fp "dcdfreq            1000"
    puts $fp "xstFreq            1000"
    puts $fp "outputEnergies     1000"
    puts $fp "outputPressure     1000"
    puts $fp ""
    puts $fp "# Protocollo di Esecuzione Matched: 1000 step Minimizzazione + 1.250.000 step (2.5 ns) NPT"
    puts $fp "minimize           1000"
    puts $fp "reinitvels         \$temperature"
    puts $fp "run                1250000"
    close $fp
}

# 1. Configura le repliche Seed 2 e Seed 3 del Wild-Type (Seed 1 e' gia' stata eseguita in 02_wt_baseline/wt_400K.dcd!)
file mkdir "04_md_simulations/WT"
write_matched_namd_conf "04_md_simulations/WT/wt_400K_seed2.conf" \
    "01_structures/2myo_wb_ionized.psf" "01_structures/2myo_wb_ionized.pdb" \
    "04_md_simulations/WT/wt_400K_seed2" 23456 71.840 58.972 54.840 7.735 -8.087 86.003
write_matched_namd_conf "04_md_simulations/WT/wt_400K_seed3.conf" \
    "01_structures/2myo_wb_ionized.psf" "01_structures/2myo_wb_ionized.pdb" \
    "04_md_simulations/WT/wt_400K_seed3" 34567 71.840 58.972 54.840 7.735 -8.087 86.003

# 2. Costruisci i sistemi CHARMM36m + TIP3P + 0.15 M NaCl per le varianti mutate
set variants_list {
    "K37L"            "03_mutations/final_models_pdb/2myo_mut_K37L.pdb"
    "T103L"           "03_mutations/final_models_pdb/2myo_mut_T103L.pdb"
    "K37L_T103L"      "03_mutations/final_models_pdb/2myo_mut_K37L_T103L.pdb"
    "H68G"            "03_mutations/final_models_pdb/2myo_mut_H68G.pdb"
    "K37L_H68G_T103L" "03_mutations/final_models_pdb/2myo_mut_K37L_H68G_T103L.pdb"
}

set summary_fp [open "04_md_simulations/matched_systems_qc_report.txt" w]
puts $summary_fp "========================================================================================================="
puts $summary_fp " REPORT CONTROLLO DI QUALITA' TOPOLOGICO (GATE 2): SISTEMI SOLVATATI MATCHED COMPARISON (CHARMM36m)"
puts $summary_fp "========================================================================================================="
puts $summary_fp [format "%-18s %-10s %-10s %-10s %-10s %-12s %-28s" "Sistema" "Atomi_Tot" "Acque" "Na+(SOD)" "Cl-(CLA)" "Carica_Tot" "Box_PBC (Lx x Ly x Lz A)"]
puts $summary_fp "---------------------------------------------------------------------------------------------------------"
puts $summary_fp [format "%-18s %-10d %-10d %-10d %-10d %-12s %-28s" "WT (2MYO)" 21461 6535 23 18 "+0.0000 e" "71.840 x 58.972 x 54.840"]

foreach {vname src_pdb} $variants_list {
    if {![file exists $src_pdb]} {
        puts "ATTENZIONE: File $src_pdb non trovato, salto $vname"
        continue
    }
    set vdir "04_md_simulations/$vname"
    file mkdir $vdir

    puts "\n=================================================================="
    puts " COSTRUZIONE SISTEMA MATCHED CHARMM36m PER VARIANTE: $vname"
    puts "=================================================================="

    mol delete all
    mol new $src_pdb

    # Assegnazione esplicita dei tautomeri verificati per le 4 Istidine (se non mutate in GLY)
    [atomselect top "resid 40"] set resname HSD
    [atomselect top "resid 67"] set resname HSE
    set r68 [atomselect top "resid 68"]
    if {[lindex [$r68 get resname] 0] in {"HIS" "H1S" "H2S" "HSD" "HSE"}} {
        $r68 set resname HSD
    }
    [atomselect top "resid 80"] set resname HSD

    # Rimozione di qualsiasi idrogeno di FoldX ("protein and noh")
    set heavy [atomselect top "protein and noh"]
    $heavy writepdb "$vdir/clean_noh.pdb"
    mol delete all

    # Generazione PSF e PDB con psfgen e topologia CHARMM36m (identico a 02_build_wt_charmm36m.tcl)
    resetpsf
    topology common/charmm36/top_all36_prot.rtf
    pdbalias atom ILE CD1 CD

    segment A {
        pdb "$vdir/clean_noh.pdb"
    }
    coordpdb "$vdir/clean_noh.pdb" A
    guesscoord

    writepdb "$vdir/${vname}_prot.pdb"
    writepsf "$vdir/${vname}_prot.psf"

    # Solvatazione TIP3P (padding 12.0 A, identico al WT)
    solvate "$vdir/${vname}_prot.psf" "$vdir/${vname}_prot.pdb" -t 12 -o "$vdir/${vname}_wb"

    # Neutralizzazione e forza ionica 0.15 M NaCl (identico al WT)
    autoionize -psf "$vdir/${vname}_wb.psf" -pdb "$vdir/${vname}_wb.pdb" -sc 0.15 -o "$vdir/${vname}_wb_ionized"

    # Misura parametri PBC e verifica carica netta = 0.0000 e
    mol delete all
    mol new "$vdir/${vname}_wb_ionized.psf"
    mol addfile "$vdir/${vname}_wb_ionized.pdb"

    set all_sys [atomselect top "all"]
    set wat_sel [atomselect top "water"]
    set sod_sel [atomselect top "resname SOD"]
    set cla_sel [atomselect top "resname CLA"]
    set tot_q   [vecsum [$all_sys get charge]]

    set minmax  [measure minmax $all_sys]
    set box_vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    set center  [measure center $all_sys]

    set lx [lindex $box_vec 0]
    set ly [lindex $box_vec 1]
    set lz [lindex $box_vec 2]
    set cx [lindex $center 0]
    set cy [lindex $center 1]
    set cz [lindex $center 2]

    set n_atoms [$all_sys num]
    set n_wat   [expr {[$wat_sel num] / 3}]
    set n_sod   [$sod_sel num]
    set n_cla   [$cla_sel num]

    puts $summary_fp [format "%-18s %-10d %-10d %-10d %-10d %+8.4f e   %.3f x %.3f x %.3f" \
        $vname $n_atoms $n_wat $n_sod $n_cla $tot_q $lx $ly $lz]

    # Genera i 3 file .conf per i 3 velocity seeds (12345, 23456, 34567)
    foreach {s_idx s_val} {1 12345 2 23456 3 34567} {
        write_matched_namd_conf "$vdir/${vname}_400K_seed${s_idx}.conf" \
            "$vdir/${vname}_wb_ionized.psf" "$vdir/${vname}_wb_ionized.pdb" \
            "$vdir/${vname}_400K_seed${s_idx}" $s_val $lx $ly $lz $cx $cy $cz
    }
    mol delete all
}

puts $summary_fp "========================================================================================================="
close $summary_fp

puts "\n=========================================================================="
puts " TUTTI I SISTEMI MATCHED CHARMM36m E I FILE .CONF MULTI-SEED SONO PRONTI!"
puts " Controlla il report di qualita' topologica (Gate 2):"
puts "   cat 04_md_simulations/matched_systems_qc_report.txt"
puts "==========================================================================\n"
exit
