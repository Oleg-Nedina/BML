# ==============================================================================
# 01_inspect_2myo.tcl
# Fase 1 (ROAD1.md - Passi 1.2 & 1.4):
# Ispezione chimico-fisica di 2MYO.pdb (sequenza, Cisteine, Istidine)
# e sovrapposizione con 3AAA_complex_ref.pdb per estrarre la Functional Mask.
#
# Esecuzione dalla cartella principale BML/:
#   - In modalita' testuale: vmd -dispdev text -e scripts/01_inspect_2myo.tcl
#   - Oppure dentro VMD (Tk Console): source scripts/01_inspect_2myo.tcl
# ==============================================================================

mol new 01_structures/raw/2MYO.pdb

puts "\n=========================================="
puts "1. SEQUENZA FASTA E COMPLETEZZA CATENA"
puts "=========================================="
set prot [atomselect top "protein"]
puts "Numero totale atomi nel PDB: [$prot num]"
set heavy [atomselect top "protein and noh"]
puts "Numero atomi pesanti (noh):  [$heavy num]"
set ca [atomselect top "protein and name CA"]
puts "Numero residui (CA):         [$ca num]"
puts "Sequenza FASTA (da VMDextensions):"
puts [getFasta $prot]

puts "\n=========================================="
puts "2. VERIFICA DELLE 3 CISTEINE (DISTANZE SG-SG)"
puts "=========================================="
set cys_sg [atomselect top "resname CYS and name SG"]
set cys_res [$cys_sg get resid]
set cys_crd [$cys_sg get {x y z}]
for {set i 0} {$i < [llength $cys_res]} {incr i} {
    for {set j [expr {$i+1}]} {$j < [llength $cys_res]} {incr j} {
        set r1 [lindex $cys_res $i]
        set r2 [lindex $cys_res $j]
        set d [vecdist [lindex $cys_crd $i] [lindex $cys_crd $j]]
        puts [format "Distanza SG(CYS %d) -- SG(CYS %d) = %.2f A" $r1 $r2 $d]
    }
}

puts "\n=========================================="
puts "3. ISPEZIONE DELLE 4 ISTIDINE (IDROGENI NMR E CONTATTI < 3.5 A)"
puts "=========================================="
foreach hres {40 67 68 80} {
    set hsel [atomselect top "resid $hres"]
    set anames [$hsel get name]
    puts "\n--- HIS $hres ---"
    puts "Atomi presenti in 2MYO.pdb per HIS $hres: $anames"
    foreach ring_n {ND1 NE2} {
        set n_sel [atomselect top "resid $hres and name $ring_n"]
        set nbrs [atomselect top "(noh and not resid $hres) and within 3.5 of (resid $hres and name $ring_n)"]
        puts "  Vicini pesanti entro 3.5 A da $ring_n:"
        foreach r [$nbrs get resname] id [$nbrs get resid] an [$nbrs get name] crd [$nbrs get {x y z}] {
            set dist [vecdist [lindex [$n_sel get {x y z}] 0] $crd]
            puts [format "    -> %s %d : %s (d = %.2f A)" $r $id $an $dist]
        }
    }
}

puts "\n=========================================="
puts "4. CONFRONTO CON 3AAA (X-RAY) E FUNCTIONAL MASK (< 4.5 A)"
puts "=========================================="
mol new 01_structures/raw/3AAA_complex_ref.pdb
# 2MYO e' mol 0, 3AAA e' mol 1
set sel_2myo [atomselect 0 "chain A and resid 2 to 118 and name CA"]
set sel_3aaa [atomselect 1 "chain C and resid 2 to 118 and name CA"]
set all_2myo [atomselect 0 "all"]
$all_2myo move [measure fit $sel_2myo $sel_3aaa]
set rmsd_val [measure rmsd $sel_2myo $sel_3aaa]
puts [format "RMSD (CA resid 2-118) tra 2MYO (Ratto NMR) e 3AAA (Uomo X-Ray): %.3f A" $rmsd_val]

set interface_sel [atomselect 1 "protein and chain C and noh and within 4.5 of (protein and chain A B and noh)"]
set mask_resids [lsort -unique -integer [$interface_sel get resid]]
puts "\nLista residui funzionali di Myotrophin a distanza < 4.5 A da Capping Protein (3AAA):"
foreach r $mask_resids {
    set r_ca [atomselect 1 "chain C and resid $r and name CA"]
    puts [format "  -> Resid %3d (%s)" $r [$r_ca get resname]]
}
puts "\nStringa di selezione per Functional Mask:"
puts "resid $mask_resids"
puts "==========================================\n"
exit
