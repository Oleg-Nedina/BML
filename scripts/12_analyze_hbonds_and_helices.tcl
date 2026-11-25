# ==============================================================================
# 12_analyze_hbonds_and_helices.tcl
# Analisi avanzata richiesta dal Prof. Gautieri (25 Nov 2026):
#   1. Rete completa di Legami a Idrogeno (-dist 3.5 -ang 30 -type unique):
#      - H-bonds totali intra-proteici (protein) -> *_hbonds35_prot.dat + details
#      - H-bonds del backbone elicoidale (protein and backbone) -> *_hbonds35_bb.dat
#   2. Matrice STRIDE residuo (1..118) x tempo (frame) -> *_stride_matrix.dat
#      (per costruire la Timeline/Heatmap della struttura secondaria)
#   3. Comportamento quantitativo di ciascuna delle 7 alpha-eliche native
#      (definite in report/01_protein_characterization.md, Yang et al. 1998):
#        H1 (3-12), H2 (15-23)          -> N-terminal (1-27)
#        H3 (39-44), H4 (49-57)         -> Ankyrin Repeat 1 (28-60)
#        H5 (71-77), H6 (82-89)         -> Ankyrin Repeat 2 (61-93)
#        H7 (112-117)                   -> C-terminal (94-118)
#      Per ogni elica esporta nel tempo:
#        - Frazione elicoidale conservata (STRIDE 'H')
#        - RMSD interno dei CA rispetto al Frame 0 (deformazione/bending)
#        - Distanza end-to-end d(CA_start - CA_end) in A (fraying/sfilacciamento!)
#
# Esecuzione dalla cartella BML/:
#   vmd -dispdev text -e scripts/12_analyze_hbonds_and_helices.tcl
# ==============================================================================

package require pbctools
package require hbonds

file mkdir "05_analysis/helices_hbonds"

set HELICES {
    {H1   3  12}
    {H2  15  23}
    {H3  39  44}
    {H4  49  57}
    {H5  71  77}
    {H6  82  89}
    {H7 112 117}
}

proc analyze_helices_hbonds {psf_file dcd_file vname seed_idx helices_list} {
    if {![file exists $dcd_file]} {
        return 0
    }
    set tag "${vname}_s${seed_idx}"
    set out_dir "05_analysis/helices_hbonds"
    set done_f "${out_dir}/${tag}_helices_geom.dat"
    if {[file exists $done_f] && ([file size $done_f] > 1000) && ([file mtime $done_f] >= [file mtime $dcd_file])} {
        puts ">>> \[SKIP CACHE\] Analisi eliche/H-bonds già aggiornata per $tag."
        return 1
    }

    puts "\n=================================================================="
    puts " ANALISI H-BONDS (3.5A/30°) E 7 ALFA-ELICHE (H1-H7): $tag"
    puts "=================================================================="

    mol delete all
    mol new $psf_file
    mol addfile $dcd_file waitfor all

    set nf [molinfo top get numframes]
    if {$nf < 5} {
        mol delete all
        return 0
    }

    # Ricomposizione PBC e allineamento sul backbone Frame 0
    catch {pbc wrap -compound residue -center com -centersel "protein" -all}
    set all [atomselect top "all"]
    set bb  [atomselect top "protein and backbone"]
    set bb0 [atomselect top "protein and backbone" frame 0]
    for {set i 0} {$i < $nf} {incr i} {
        $all frame $i
        $bb frame $i
        $all move [measure fit $bb $bb0]
    }

    # 1. H-Bonds completi (criterio standard letteratura: dist 3.5 A, angolo 30 deg, -type unique)
    set prot [atomselect top "protein"]
    file delete -force "${out_dir}/${tag}_hbonds35_prot.dat" "${out_dir}/${tag}_hbonds35_prot_details.txt"
    hbonds -sel1 $prot -upsel no -writefile yes -polar yes -plot no \
           -dist 3.5 -ang 30 -type unique -frames 0:1:[expr {$nf - 1}] \
           -outfile "${out_dir}/${tag}_hbonds35_prot.dat" \
           -detailout "${out_dir}/${tag}_hbonds35_prot_details.txt"

    file delete -force "${out_dir}/${tag}_hbonds35_bb.dat"
    hbonds -sel1 $bb -upsel no -writefile yes -polar yes -plot no \
           -dist 3.5 -ang 30 -type unique -frames 0:1:[expr {$nf - 1}] \
           -outfile "${out_dir}/${tag}_hbonds35_bb.dat"

    # 2. Prepara selezioni per le 7 eliche native (CA, estremità N/C)
    set ca_all [atomselect top "protein and name CA"]
    foreach h $helices_list {
        lassign $h hname r1 r2
        set h_ca($hname)   [atomselect top "protein and name CA and resid $r1 to $r2"]
        set h_ca0($hname)  [atomselect top "protein and name CA and resid $r1 to $r2" frame 0]
        set h_n($hname)    [atomselect top "protein and name CA and resid $r1"]
        set h_c($hname)    [atomselect top "protein and name CA and resid $r2"]
        set h_len($hname)  [expr {$r2 - $r1 + 1}]
    }

    set fp_mat  [open "${out_dir}/${tag}_stride_matrix.dat" w]
    set fp_geom [open "${out_dir}/${tag}_helices_geom.dat" w]
    puts $fp_geom "# Frame H1_frac H1_rmsd H1_end2end H2_frac H2_rmsd H2_end2end H3_frac H3_rmsd H3_end2end H4_frac H4_rmsd H4_end2end H5_frac H5_rmsd H5_end2end H6_frac H6_rmsd H6_end2end H7_frac H7_rmsd H7_end2end"

    for {set f 0} {$f < $nf} {incr f} {
        animate goto $f
        mol ssrecalc top
        $ca_all frame $f
        set ss_list [$ca_all get structure]

        # Salva riga matrice STRIDE (1 = H/G/I elica, 0 = altro) ogni 2 frame per leggerezza
        if {$f % 2 == 0} {
            set num_row {}
            foreach c $ss_list {
                if {$c eq "H"} {
                    lappend num_row 2
                } elseif {$c eq "G" || $c eq "I"} {
                    lappend num_row 1
                } else {
                    lappend num_row 0
                }
            }
            puts $fp_mat "$f [join $num_row " "]"
        }

        # Calcola metriche per ciascuna delle 7 eliche H1..H7
        set line_out [format "%d" $f]
        foreach h $helices_list {
            lassign $h hname r1 r2
            $h_ca($hname) frame $f
            $h_n($hname) frame $f
            $h_c($hname) frame $f

            # Frazione di residui in stato alpha-helix ('H')
            set sub_ss [$h_ca($hname) get structure]
            set n_h [regexp -all -- H $sub_ss]
            set frac_h [expr {1.0 * $n_h / $h_len($hname)}]

            # RMSD interno della singola elica (dopo fit locale dell'elica su se stessa al frame 0)
            set tmat [measure fit $h_ca($hname) $h_ca0($hname)]
            set coords_orig [$h_ca($hname) get {x y z}]
            $h_ca($hname) move $tmat
            set h_rmsd [measure rmsd $h_ca($hname) $h_ca0($hname)]
            $h_ca($hname) set {x y z} $coords_orig

            # Distanza end-to-end (CA_r1 -- CA_r2) in A
            set p1 [lindex [$h_n($hname) get {x y z}] 0]
            set p2 [lindex [$h_c($hname) get {x y z}] 0]
            set d_ee [vecdist $p1 $p2]

            append line_out [format " %.4f %.4f %.4f" $frac_h $h_rmsd $d_ee]
        }
        puts $fp_geom $line_out
    }

    close $fp_mat
    close $fp_geom
    mol delete all
    return 1
}

# 1. Wild-Type (Seed 1, 2, 3 — preferisce wt_400K_seed1.dcd con seed 12345 se presente)
set wt_s1_dcd "04_md_simulations/WT/wt_400K_seed1.dcd"
if {![file exists $wt_s1_dcd]} { set wt_s1_dcd "02_wt_baseline/wt_400K.dcd" }
analyze_helices_hbonds "01_structures/2myo_wb_ionized.psf" $wt_s1_dcd "WT" 1 $HELICES
analyze_helices_hbonds "01_structures/2myo_wb_ionized.psf" "04_md_simulations/WT/wt_400K_seed2.dcd" "WT" 2 $HELICES
analyze_helices_hbonds "01_structures/2myo_wb_ionized.psf" "04_md_simulations/WT/wt_400K_seed3.dcd" "WT" 3 $HELICES

# 2. Varianti Mutate (Seed 1, 2, 3)
foreach vname {"K37L" "T103L" "K37L_T103L" "K37L_H68G_T103L" "H68G"} {
    set psf_f "04_md_simulations/${vname}/${vname}_wb_ionized.psf"
    if {![file exists $psf_f]} { continue }
    for {set s 1} {$s <= 3} {incr s} {
        set dcd_f "04_md_simulations/${vname}/${vname}_400K_seed${s}.dcd"
        analyze_helices_hbonds $psf_f $dcd_f $vname $s $HELICES
    }
}

puts "\n=========================================================================="
puts " ANALISI H-BONDS (3.5A/30°) E 7 ALFA-ELICHE COMPLETATA!"
puts " Esegui ora: python3 scripts/12_plot_hbonds_and_helices.py"
puts "==========================================================================\n"
exit
