# ==============================================================================
# 03_analyze_wt_baseline.tcl
# Fase 2 (ROAD2.md - Passo 2.4):
# Integra al 100% lo script ufficiale del Prof. Gautieri (laboratorio/analysisMD.tcl)
# con l'analisi dei Contatti Nativi (Q-value, Bekker et al. 2019):
#
#   [Da analysisMD.tcl del Prof. Gautieri]:
#   1. Allineamento della traiettoria sul backbone del Frame 0 (measure fit)
#   2. RMSD del backbone nel tempo (rmsd_${tag}.dat)
#   3. RMSF dei Carboni Alfa (CA) per residuo (rmsfCA_${tag}.dat)
#   4. Legami a Idrogeno (hbonds_${tag}.dat + dettaglio occupancy hbonds_details_${tag}.txt)
#   5. SASA - Solvent Accessible Surface Area con raggio sonda 1.4 A (sasa_${tag}.dat)
#   6. Rgyr - Raggio di Girazione della proteina (rgyr_${tag}.dat)
#   7. Alpha-Helicity (%) tramite STRIDE / mol ssrecalc (alpha_${tag}.dat)
#
#   [Da Tutorials/nativecontacts/native_contacts.tcl & Bekker et al., 2019]:
#   8. Q-value nel tempo: all-all, hydrophilic-all, hydrophobic-small e ANK1..ANK4 (q_time_${tag}.dat)
#   9. <Q_i> per singolo residuo (1..118) sugli ultimi 30% dei frame a 400 K,
#      report residui instabili (<Q_i> < 0.60) e PDB colorato per Beta (wt_400K_q_beta.pdb).
#
# Esecuzione dalla cartella principale BML/:
#   vmd -dispdev text -e scripts/03_analyze_wt_baseline.tcl
# ==============================================================================

source scripts/VMDextensions.tcl
package require pbctools
package require hbonds

# Rileva automaticamente se stiamo validando (smoke test) o analizzando la run notturna
if {[file exists "02_wt_baseline/wt_400K.dcd"] && [file exists "02_wt_baseline/wt_300K.dcd"]} {
    set dcd_300 "02_wt_baseline/wt_300K.dcd"
    set dcd_400 "02_wt_baseline/wt_400K.dcd"
    puts ">>> MODALITA' PRODUZIONE: Analisi di wt_300K.dcd e wt_400K.dcd"
} else {
    set dcd_300 "02_wt_baseline/benchmark_out.dcd"
    set dcd_400 "02_wt_baseline/smoke_400K.dcd"
    puts ">>> MODALITA' SMOKE-TEST VALIDAZIONE: Analisi di benchmark_out.dcd e smoke_400K.dcd"
}

# Lista dei 24 residui della Functional Mask (verificata al Passo 1.4 da 3AAA)
set functional_mask {4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100}

# Procedura completa che unisce analysisMD.tcl + Q-value per una traiettoria (300K o 400K)
proc analyze_traj {psf_file pdb_file dcd_file tag} {
    mol new $psf_file
    mol addfile $pdb_file
    mol addfile $dcd_file waitfor all

    set nf [molinfo top get numframes]
    set step 1
    puts "\n=================================================================="
    puts " ANALISI TRAIETTORIA $tag ($dcd_file): $nf frames totali"
    puts "=================================================================="

    # Selezioni standard da analysisMD.tcl del Prof. Gautieri
    set all       [atomselect top all]
    set bb        [atomselect top "protein and backbone"]
    set bb0       [atomselect top "protein and backbone" frame 0]
    set protein   [atomselect top "protein"]
    set ca        [atomselect top "protein and name CA"]
    set numres    [$ca num]

    # 1. Allineamento della traiettoria sul backbone del frame 0 (analysisMD.tcl righe 17-27)
    puts " -> (1/7) Allineamento traiettoria sul backbone (Frame 0)..."
    set reference [atomselect top "protein and backbone" frame 0]
    set compare   [atomselect top "protein and backbone"]
    for {set i 0} {$i < $nf} {incr i} {
        $all frame $i
        $compare frame $i
        set trans_mat [measure fit $compare $reference]
        $all move $trans_mat
    }

    # 2. Calcolo RMSD del backbone (analysisMD.tcl righe 29-41)
    puts " -> (2/7) Calcolo RMSD sul backbone..."
    set out_rmsd [open "02_wt_baseline/rmsd_${tag}.dat" w]
    puts $out_rmsd "# Frame RMSD_Backbone(A)"
    for {set i 0} {$i < $nf} {incr i $step} {
        $bb frame $i
        $bb0 frame 0
        set rmsd [measure rmsd $bb0 $bb]
        puts $out_rmsd [format "%d %.4f" $i $rmsd]
    }
    close $out_rmsd

    # 3. Calcolo RMSF sui CA (analysisMD.tcl righe 44-54)
    puts " -> (3/7) Calcolo RMSF sui CA..."
    set start_f [expr {$nf > 2 ? 1 : 0}]
    set rmsfCA [measure rmsf $ca first $start_f last [expr {$nf - 1}]]
    set firstAA [lindex [$ca get resid] 0]
    set out_rmsf [open "02_wt_baseline/rmsf_${tag}.dat" w]
    puts $out_rmsf "# Resid RMSF_CA(A)"
    for {set i 0} {$i < [$ca num]} {incr i} {
        puts $out_rmsf [format "%d %.4f" [expr {$i + $firstAA}] [lindex $rmsfCA $i]]
    }
    close $out_rmsf

    # 4. Calcolo Hydrogen Bonds (analysisMD.tcl righe 57-62 + dettaglio coppie donatore-accettore)
    puts " -> (4/7) Calcolo Hydrogen Bonds (backbone e proteina completa)..."
    file delete -force "02_wt_baseline/hbonds_${tag}.dat" "02_wt_baseline/hbonds_details_${tag}.txt"
    hbonds -sel1 $protein -upsel no -writefile yes -polar yes -plot no \
           -frames 0:$step:[expr {$nf - 1}] \
           -outfile "02_wt_baseline/hbonds_${tag}.dat" \
           -detailout "02_wt_baseline/hbonds_details_${tag}.txt"

    # 5. Calcolo SASA e Raggio di Girazione Rgyr (analysisMD.tcl righe 65-90)
    puts " -> (5/7) Calcolo SASA (r_probe = 1.4 A) e Raggio di Girazione (Rgyr)..."
    set out_sasa [open "02_wt_baseline/sasa_${tag}.dat" w]
    set out_rgyr [open "02_wt_baseline/rgyr_${tag}.dat" w]
    puts $out_sasa "# Frame SASA(A^2)"
    puts $out_rgyr "# Frame Rgyr(A)"
    for {set i 0} {$i < $nf} {incr i $step} {
        $protein frame $i
        set sasa [measure sasa 1.4 $protein]
        set rgyr [measure rgyr $protein]
        puts $out_sasa [format "%d %.2f" $i $sasa]
        puts $out_rgyr [format "%d %.4f" $i $rgyr]
    }
    close $out_sasa
    close $out_rgyr

    # 6. Calcolo Frazione di Alpha-Helix e Beta-Sheet (analysisMD.tcl righe 92-125)
    puts " -> (6/7) Calcolo contenuto di Struttura Secondaria (Alpha-Helix %)..."
    set out_alpha [open "02_wt_baseline/alpha_${tag}.dat" w]
    puts $out_alpha "# Frame AlphaPerc AlphaResidues"
    for {set i 0} {$i < $nf} {incr i $step} {
        animate goto $i
        mol ssrecalc top
        $ca frame $i
        set st [$ca get structure]
        set Alpha [regexp -all -- H $st]
        set AlphaPerc [expr {1.00 * $Alpha / $numres}]
        puts $out_alpha [format "%d %.4f %d" $i $AlphaPerc $Alpha]
    }
    close $out_alpha

    # 7. Calcolo Q-value nel tempo (Globale, Chimico Bekker et al. e per Ankyrin Repeat ANK1..4)
    puts " -> (7/7) Calcolo Q-value (all-all, hydrophilic-all, hydrophobic-small, ANK1..4)..."
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

    set n_all  [llength $ref_all]
    set n_phil [llength $ref_phil]
    set n_ps   [llength $ref_ps]
    set n_a1   [llength $ref_ank1]
    set n_a2   [llength $ref_ank2]
    set n_a3   [llength $ref_ank3]
    set n_a4   [llength $ref_ank4]

    set out_q [open "02_wt_baseline/q_time_${tag}.dat" w]
    puts $out_q "# Frame Q_all Q_phil_all Q_phob_small Q_ANK1 Q_ANK2 Q_ANK3 Q_ANK4"

    for {set f 0} {$f < $nf} {incr f} {
        foreach s [list $sel_all $sel_phil $sel_phob $sel_small $sel_ank1 $sel_ank2 $sel_ank3 $sel_ank4] {
            $s frame $f
        }
        set q_all  [expr {1.0 * [measureNativeContacts $ref_all  4.5 $sel_all] / $n_all}]
        set q_phil [expr {1.0 * [measureNativeContacts $ref_phil 4.5 $sel_phil $sel_all] / $n_phil}]
        set q_ps   [expr {1.0 * [measureNativeContacts $ref_ps   4.5 $sel_phob $sel_small] / $n_ps}]
        set q_a1   [expr {1.0 * [measureNativeContacts $ref_ank1 4.5 $sel_ank1 $sel_all] / $n_a1}]
        set q_a2   [expr {1.0 * [measureNativeContacts $ref_ank2 4.5 $sel_ank2 $sel_all] / $n_a2}]
        set q_a3   [expr {1.0 * [measureNativeContacts $ref_ank3 4.5 $sel_ank3 $sel_all] / $n_a3}]
        set q_a4   [expr {1.0 * [measureNativeContacts $ref_ank4 4.5 $sel_ank4 $sel_all] / $n_a4}]

        puts $out_q [format "%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f" $f $q_all $q_phil $q_ps $q_a1 $q_a2 $q_a3 $q_a4]
    }
    close $out_q
    mol delete top
}

# Esegui analisi completa su 300K e 400K
analyze_traj "01_structures/2myo_wb_ionized.psf" "01_structures/2myo_wb_ionized.pdb" $dcd_300 "300K"
analyze_traj "01_structures/2myo_wb_ionized.psf" "01_structures/2myo_wb_ionized.pdb" $dcd_400 "400K"

# ==============================================================================
# 8. Calcolo <Q_i> per singolo residuo (1..118) negli ultimi 30% dei frame a 400K
#    (Metodo Bekker, Ma & Kamiya 2019, Protein Science 28:429-438, Table S1)
# ==============================================================================
mol new "01_structures/2myo_wb_ionized.psf"
mol addfile "01_structures/2myo_wb_ionized.pdb"
mol addfile $dcd_400 waitfor all

set nf [molinfo top get numframes]
set f_start [expr {int($nf * 0.70)}]
if {$f_start < 1} { set f_start 1 }
set n_eval [expr {$nf - $f_start}]

# Allinea prima di salvare il PDB finale
set all_at [atomselect top all]
set ref_bb [atomselect top "protein and backbone" frame 0]
set cur_bb [atomselect top "protein and backbone"]
for {set f 0} {$f < $nf} {incr f} {
    $all_at frame $f
    $cur_bb frame $f
    $all_at move [measure fit $cur_bb $ref_bb]
}

puts "\n=================================================================="
puts " CALCOLO <Q_i> PER RESIDUO A 400K (Frame $f_start .. [expr {$nf-1}])"
puts "=================================================================="

set out_qi [open "02_wt_baseline/q_per_residue_400K.dat" w]
puts $out_qi "# Resid Resname N_contacts_0 Mean_Qi Is_Functional_Mask"

set out_rep [open "02_wt_baseline/unstable_residues_report.txt" w]
puts $out_rep "=========================================================================="
puts $out_rep " REPORT RESIDUI LOCALMENTE INSTABILI A 400 K (Bekker et al., 2019)"
puts $out_rep " Finestra temporale analizzata: Frame $f_start - [expr {$nf-1}] (su $nf)"
puts $out_rep "=========================================================================="
puts $out_rep [format "%-6s %-8s %-8s %-10s %-12s %-20s" "Resid" "Resname" "Repeat" "N_cont(0)" "<Q_i>_400K" "Stato / Candidato"]
puts $out_rep "--------------------------------------------------------------------------"

set prot_all [atomselect top "protein"]
$prot_all set beta 1.0

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

    if {$r <= 34} { set rep "ANK1" } elseif {$r <= 67} { set rep "ANK2" } elseif {$r <= 100} { set rep "ANK3" } else { set rep "ANK4" }

    set in_mask [expr {[lsearch -exact $functional_mask $r] >= 0 ? 1 : 0}]
    if {$in_mask} {
        set status "MASKED (CP-binding)"
    } elseif {$r <= 3 || $r >= 116} {
        set status "TERMINAL"
    } elseif {$mean_qi < 0.60 && $nnc_r >= 30} {
        set status ">>> UNSTABLE (Q<0.6) <<<"
    } elseif {$mean_qi < 0.72} {
        set status "WEAK (0.6<=Q<0.72)"
    } else {
        set status "STABLE"
    }

    puts $out_qi [format "%d %s %d %.4f %d" $r $rname $nnc_r $mean_qi $in_mask]
    puts $out_rep [format "%-6d %-8s %-8s %-10d %-12.4f %-20s" $r $rname $rep $nnc_r $mean_qi $status]

    set r_all_atoms [atomselect top "protein and resid $r" frame [expr {$nf - 1}]]
    $r_all_atoms set beta $mean_qi

    $r_sel delete; $env_sel delete; $r_ca delete; $r_all_atoms delete
}

close $out_qi
close $out_rep

set last_prot [atomselect top "protein" frame [expr {$nf - 1}]]
$last_prot writepdb "02_wt_baseline/wt_400K_q_beta.pdb"

set init_prot [atomselect top "protein" frame 0]
$init_prot writepdb "02_wt_baseline/wt_initial_ref.pdb"

puts "=================================================================="
puts " Analisi completata! Tutti i file (analysisMD.tcl + Q-value) sono"
puts " stati salvati in 02_wt_baseline/:"
puts "   - rmsd_300K.dat, rmsd_400K.dat"
puts "   - rmsf_300K.dat, rmsf_400K.dat"
puts "   - hbonds_300K.dat, hbonds_400K.dat (+ hbonds_details_*.txt)"
puts "   - sasa_300K.dat, sasa_400K.dat"
puts "   - rgyr_300K.dat, rgyr_400K.dat"
puts "   - alpha_300K.dat, alpha_400K.dat"
puts "   - q_time_300K.dat, q_time_400K.dat"
puts "   - q_per_residue_400K.dat & unstable_residues_report.txt"
puts "==================================================================\n"
exit
