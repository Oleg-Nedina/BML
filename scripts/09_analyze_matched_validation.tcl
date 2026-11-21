# ==============================================================================
# 09_analyze_matched_validation.tcl
# Punto 5 della Scaletta del Prof. Gautieri (Laboratory_Project_2026.pdf, Slide 15, 22, 23)
# & Livelli 3, 4, 5 della Politica di Validazione (ROAD4_ROAD5_POLITICA_DI_VALIDAZIONE.md):
#
# Analizza tutte le repliche multi-seed a 400 K disponibili per:
#   - WT         (Seed 1 in 02_wt_baseline/wt_400K.dcd + Seed 2, 3 in 04_md_simulations/WT/)
#   - K37L       (Seed 1, 2, 3 in 04_md_simulations/K37L/)
#   - T103L      (Seed 1, 2, 3 in 04_md_simulations/T103L/)
#   - K37L_T103L (Seed 1, 2, 3 in 04_md_simulations/K37L_T103L/)
#   - (Opzionali: H68G e K37L_H68G_T103L se simulate)
#
# Per ogni replica calcola ed esporta in 05_analysis/<variant>_seed<k>_*:
#   1. Ricomposizione PBC (pbc wrap) e allineamento sul backbone Frame 0 (Slide 22)
#   2. Backbone RMSD(t) e Functional Mask RMSD(t) sui 24 residui di legame 3AAA
#   3. RMSF(CA) per residuo (1..118)
#   4. Hydrogen Bonds totali nel tempo + dettaglio occupancy (-type unique)
#   5. SASA(t), Raggio di Girazione Rgyr(t) e % Alpha-Helix STRIDE (analysisMD.tcl)
#   6. Q-value nel tempo (Bekker et al. 2019 Table II + Lowe & Itzhaki 2007):
#      - Q_all, Q_phil_all (r=0.84 con Tm!), Q_phob_small
#      - Q_ANK1, Q_ANK2, Q_ANK3, Q_ANK4
#      - Q_int12 (interfaccia ANK1-ANK2), Q_int23 (ANK2-ANK3), Q_int34 (ANK3-ANK4)
#   7. <Q_i>_400K per singolo residuo (1..118) nella finestra equilibrata finale (ultimi 30% frame)
#      e conteggio dei residui dello scaffold localmente instabili (<Q_i> < 0.60, N_cont >= 30).
#
# Esecuzione dalla cartella principale BML/:
#   vmd -dispdev text -e scripts/09_analyze_matched_validation.tcl
# ==============================================================================

source scripts/VMDextensions.tcl
package require pbctools
package require hbonds

file mkdir "05_analysis"

set functional_mask {4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100}

proc analyze_replica {psf_file pdb_file dcd_file vname seed_idx functional_mask} {
    if {![file exists $dcd_file]} {
        return 0
    }
    set tag "${vname}_s${seed_idx}"
    set done_f "05_analysis/${tag}_q_per_residue.dat"
    if {[file exists $done_f] && ([file size $done_f] > 500) && ([file mtime $done_f] >= [file mtime $dcd_file])} {
        puts ">>> \[SKIP CACHE\] Replica $tag già analizzata in 05_analysis/ ($done_f aggiornato)."
        return 1
    }
    puts "\n=================================================================="
    puts " ANALISI REPLICA MATCHED: $tag ($dcd_file)"
    puts "=================================================================="

    mol delete all
    mol new $psf_file
    mol addfile $pdb_file
    mol addfile $dcd_file waitfor all

    set nf [molinfo top get numframes]
    if {$nf < 5} {
        puts "ATTENZIONE: Traiettoria $dcd_file ha solo $nf frame, salto."
        mol delete all
        return 0
    }

    set all     [atomselect top all]
    set bb      [atomselect top "protein and backbone"]
    set bb0     [atomselect top "protein and backbone" frame 0]
    set mask_bb [atomselect top "protein and backbone and resid $functional_mask"]
    set mask_0  [atomselect top "protein and backbone and resid $functional_mask" frame 0]
    set protein [atomselect top "protein"]
    set ca      [atomselect top "protein and name CA"]
    set numres  [$ca num]

    # 1. Ricomposizione PBC e allineamento sul backbone Frame 0 (Slide 22)
    catch {pbc wrap -compound residue -center com -centersel "protein" -all}
    for {set i 0} {$i < $nf} {incr i} {
        $all frame $i
        $bb frame $i
        $all move [measure fit $bb $bb0]
    }

    # 2. RMSD Backbone globale + RMSD Functional Mask (24 residui 3AAA), SASA, Rgyr, % Alpha-Helix
    set out_obs [open "05_analysis/${tag}_observables.dat" w]
    puts $out_obs "# Frame RMSD_BB(A) RMSD_Mask3AAA(A) Rgyr(A) SASA(A^2) AlphaPerc AlphaRes"
    for {set i 0} {$i < $nf} {incr i} {
        animate goto $i
        mol ssrecalc top
        $bb frame $i
        $mask_bb frame $i
        $protein frame $i
        $ca frame $i

        set rmsd_bb   [measure rmsd $bb0 $bb]
        set rmsd_mask [measure rmsd $mask_0 $mask_bb]
        set rgyr_val  [measure rgyr $protein]
        set sasa_val  [measure sasa 1.4 $protein]
        set st        [$ca get structure]
        set n_alpha   [regexp -all -- H $st]
        set p_alpha   [expr {1.00 * $n_alpha / $numres}]

        puts $out_obs [format "%d %.4f %.4f %.4f %.2f %.4f %d" \
            $i $rmsd_bb $rmsd_mask $rgyr_val $sasa_val $p_alpha $n_alpha]
    }
    close $out_obs

    # 3. RMSF sui Carboni Alfa (CA)
    set start_f [expr {$nf > 2 ? 1 : 0}]
    set rmsfCA [measure rmsf $ca first $start_f last [expr {$nf - 1}]]
    set firstAA [lindex [$ca get resid] 0]
    set out_rmsf [open "05_analysis/${tag}_rmsf.dat" w]
    puts $out_rmsf "# Resid RMSF_CA(A)"
    for {set i 0} {$i < $numres} {incr i} {
        puts $out_rmsf [format "%d %.4f" [expr {$i + $firstAA}] [lindex $rmsfCA $i]]
    }
    close $out_rmsf

    # 4. Hydrogen Bonds (-type unique)
    file delete -force "05_analysis/${tag}_hbonds.dat" "05_analysis/${tag}_hbonds_details.txt"
    hbonds -sel1 $protein -upsel no -writefile yes -polar yes -plot no \
           -type unique -frames 0:1:[expr {$nf - 1}] \
           -outfile "05_analysis/${tag}_hbonds.dat" \
           -detailout "05_analysis/${tag}_hbonds_details.txt"

    # 5. Q-value nel tempo (Classi Chimiche Bekker + Singoli Repeat ANK1..4 + Interfacce Inter-Repeat)
    set sel_all   [atomselect top "protein and noh"]
    set sel_phil  [atomselect top "protein and noh and resname ASP GLU GLN ASN ARG LYS HIS HSD HSE HSP"]
    set sel_phob  [atomselect top "protein and noh and resname PHE TYR TRP LEU VAL ILE MET CYS PRO"]
    set sel_small [atomselect top "protein and noh and resname GLY ALA SER THR"]
    set sel_ank1  [atomselect top "protein and noh and resid 1 to 34"]
    set sel_ank2  [atomselect top "protein and noh and resid 35 to 67"]
    set sel_ank3  [atomselect top "protein and noh and resid 68 to 100"]
    set sel_ank4  [atomselect top "protein and noh and resid 101 to 118"]

    foreach s [list $sel_all $sel_phil $sel_phob $sel_small $sel_ank1 $sel_ank2 $sel_ank3 $sel_ank4] {
        $s frame 0
    }

    set ref_all   [prepareNativeContacts 4.5 $sel_all]
    set ref_phil  [prepareNativeContacts 4.5 $sel_phil $sel_all]
    set ref_ps    [prepareNativeContacts 4.5 $sel_phob $sel_small]
    set ref_ank1  [prepareNativeContacts 4.5 $sel_ank1 $sel_all]
    set ref_ank2  [prepareNativeContacts 4.5 $sel_ank2 $sel_all]
    set ref_ank3  [prepareNativeContacts 4.5 $sel_ank3 $sel_all]
    set ref_ank4  [prepareNativeContacts 4.5 $sel_ank4 $sel_all]
    set ref_int12 [prepareNativeContacts 4.5 $sel_ank1 $sel_ank2]
    set ref_int23 [prepareNativeContacts 4.5 $sel_ank2 $sel_ank3]
    set ref_int34 [prepareNativeContacts 4.5 $sel_ank3 $sel_ank4]

    set n_all   [llength $ref_all]
    set n_phil  [llength $ref_phil]
    set n_ps    [llength $ref_ps]
    set n_a1    [llength $ref_ank1]
    set n_a2    [llength $ref_ank2]
    set n_a3    [llength $ref_ank3]
    set n_a4    [llength $ref_ank4]
    set n_i12   [expr {[llength $ref_int12] > 0 ? [llength $ref_int12] : 1}]
    set n_i23   [expr {[llength $ref_int23] > 0 ? [llength $ref_int23] : 1}]
    set n_i34   [expr {[llength $ref_int34] > 0 ? [llength $ref_int34] : 1}]

    set out_q [open "05_analysis/${tag}_q_time.dat" w]
    puts $out_q "# Frame Q_all Q_phil_all Q_phob_small Q_ANK1 Q_ANK2 Q_ANK3 Q_ANK4 Q_int12 Q_int23 Q_int34"
    for {set f 0} {$f < $nf} {incr f} {
        foreach s [list $sel_all $sel_phil $sel_phob $sel_small $sel_ank1 $sel_ank2 $sel_ank3 $sel_ank4] {
            $s frame $f
        }
        set q_all  [expr {1.0 * [measureNativeContacts $ref_all   4.5 $sel_all] / $n_all}]
        set q_phil [expr {1.0 * [measureNativeContacts $ref_phil  4.5 $sel_phil $sel_all] / $n_phil}]
        set q_ps   [expr {1.0 * [measureNativeContacts $ref_ps    4.5 $sel_phob $sel_small] / $n_ps}]
        set q_a1   [expr {1.0 * [measureNativeContacts $ref_ank1  4.5 $sel_ank1 $sel_all] / $n_a1}]
        set q_a2   [expr {1.0 * [measureNativeContacts $ref_ank2  4.5 $sel_ank2 $sel_all] / $n_a2}]
        set q_a3   [expr {1.0 * [measureNativeContacts $ref_ank3  4.5 $sel_ank3 $sel_all] / $n_a3}]
        set q_a4   [expr {1.0 * [measureNativeContacts $ref_ank4  4.5 $sel_ank4 $sel_all] / $n_a4}]
        set q_i12  [expr {1.0 * [measureNativeContacts $ref_int12 4.5 $sel_ank1 $sel_ank2] / $n_i12}]
        set q_i23  [expr {1.0 * [measureNativeContacts $ref_int23 4.5 $sel_ank2 $sel_ank3] / $n_i23}]
        set q_i34  [expr {1.0 * [measureNativeContacts $ref_int34 4.5 $sel_ank3 $sel_ank4] / $n_i34}]

        puts $out_q [format "%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f" \
            $f $q_all $q_phil $q_ps $q_a1 $q_a2 $q_a3 $q_a4 $q_i12 $q_i23 $q_i34]
    }
    close $out_q

    # 6. <Q_i>_400K per singolo residuo sugli ultimi 30% dei frame (Bekker et al. 2019)
    set f_start [expr {int($nf * 0.70)}]
    if {$f_start < 1} { set f_start 1 }
    set n_eval [expr {$nf - $f_start}]

    set out_qi [open "05_analysis/${tag}_q_per_residue.dat" w]
    puts $out_qi "# Resid Resname N_contacts_0 Mean_Qi Is_Functional_Mask"

    for {set r 1} {$r <= 118} {incr r} {
        set r_sel   [atomselect top "protein and noh and resid $r"]
        set env_sel [atomselect top "protein and noh and not resid $r"]
        $r_sel frame 0
        $env_sel frame 0

        set r_ca [atomselect top "protein and resid $r and name CA" frame 0]
        set rname [$r_ca get resname]

        set ref_r [prepareNativeContacts 4.5 $r_sel $env_sel]
        set nnc_r [llength $ref_r]

        if {$nnc_r > 0} {
            set sum_q 0.0
            for {set f $f_start} {$f < $nf} {incr f} {
                $r_sel frame $f
                $env_sel frame $f
                set nc [measureNativeContacts $ref_r 4.5 $r_sel $env_sel]
                set sum_q [expr {$sum_q + (1.0 * $nc / $nnc_r)}]
            }
            set mean_qi [expr {$sum_q / $n_eval}]
        } else {
            set mean_qi 1.0
        }

        set in_mask [expr {[lsearch -exact $functional_mask $r] >= 0 ? 1 : 0}]
        puts $out_qi [format "%d %s %d %.4f %d" $r $rname $nnc_r $mean_qi $in_mask]

        $r_sel delete; $env_sel delete; $r_ca delete
    }
    close $out_qi

    mol delete all
    return 1
}

# 1. Analizza le 3 repliche del Wild-Type (Seed 1 in 02_wt_baseline/wt_400K.dcd, Seed 2 e 3 in 04_md_simulations/WT/)
analyze_replica "01_structures/2myo_wb_ionized.psf" "01_structures/2myo_wb_ionized.pdb" \
    "02_wt_baseline/wt_400K.dcd" "WT" 1 $functional_mask
analyze_replica "01_structures/2myo_wb_ionized.psf" "01_structures/2myo_wb_ionized.pdb" \
    "04_md_simulations/WT/wt_400K_seed2.dcd" "WT" 2 $functional_mask
analyze_replica "01_structures/2myo_wb_ionized.psf" "01_structures/2myo_wb_ionized.pdb" \
    "04_md_simulations/WT/wt_400K_seed3.dcd" "WT" 3 $functional_mask

# 2. Analizza le repliche delle varianti mutate in 04_md_simulations/<vname>/
foreach vname {"K37L" "T103L" "K37L_T103L" "H68G" "K37L_H68G_T103L"} {
    set psf_f "04_md_simulations/${vname}/${vname}_wb_ionized.psf"
    set pdb_f "04_md_simulations/${vname}/${vname}_wb_ionized.pdb"
    if {![file exists $psf_f]} { continue }
    for {set s 1} {$s <= 3} {incr s} {
        set dcd_f "04_md_simulations/${vname}/${vname}_400K_seed${s}.dcd"
        analyze_replica $psf_f $pdb_f $dcd_f $vname $s $functional_mask
    }
}

puts "\n=========================================================================="
puts " ANALISI TRAIETTORIE MATCHED COMPLETATA! File salvati in 05_analysis/"
puts " Esegui ora: python3 scripts/10_plot_final_validation_bekker.py"
puts "==========================================================================\n"
exit
