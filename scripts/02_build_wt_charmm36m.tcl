# ==============================================================================
# 02_build_wt_charmm36m.tcl
# Fase 1 (ROAD1.md - Passi 1.5 & 1.6):
# 1. Rimozione idrogeni NMR ("protein and noh") da 2MYO.pdb
# 2. Assegnazione razionale dei tautomeri delle 4 Istidine (verificati al Passo 1.2):
#    - HIS 40 -> HSD (ND1-HD1 dona H-bond a SER 74 OG, d = 2.44 A)
#    - HIS 67 -> HSE (ND1 accetta da LYS 37 NZ, d = 3.24 A; NE2-HE2 dona a GLY 35 O, d = 2.40 A)
#    - HIS 68 -> HSD (esposta al solvente, tautomero NMR nativo)
#    - HIS 80 -> HSD (contatto apolare con LEU 48, tautomero NMR nativo)
# 3. Generazione topologia CHARMM36m (top_all36_prot.rtf) -> 2myo_wt.psf / .pdb
# 4. Solvatazione esplicita TIP3P (padding 12 A) -> 2myo_wb.psf / .pdb
# 5. Neutralizzazione e forza ionica fisiologica (0.15 M NaCl) -> 2myo_wb_ionized.psf / .pdb
# 6. Calcolo esatto di cellBasisVector1/2/3 e cellOrigin per NAMD.
#
# Esecuzione dalla cartella principale BML/:
#   vmd -dispdev text -e scripts/02_build_wt_charmm36m.tcl
# ==============================================================================

mol new 01_structures/raw/2MYO.pdb

# 1. Assegnazione esplicita dei tautomeri verificati per le 4 Istidine
[atomselect top "resid 40"] set resname HSD
[atomselect top "resid 67"] set resname HSE
[atomselect top "resid 68"] set resname HSD
[atomselect top "resid 80"] set resname HSD

# 2. Estrazione dei soli atomi pesanti della proteina (senza idrogeni NMR)
set heavy [atomselect top "protein and noh"]
$heavy writepdb 01_structures/2myo_clean_noh.pdb
mol delete top

# 3. Generazione di PSF e PDB del Wild-Type con psfgen e CHARMM36m
package require psfgen
resetpsf
topology common/charmm36/top_all36_prot.rtf
pdbalias atom ILE CD1 CD

segment A {
    pdb 01_structures/2myo_clean_noh.pdb
}
coordpdb 01_structures/2myo_clean_noh.pdb A
guesscoord

writepdb 01_structures/2myo_wt.pdb
writepsf 01_structures/2myo_wt.psf

# Controllo di qualita' su 2myo_wt (verifica carica e coordinate)
mol new 01_structures/2myo_wt.psf
mol addfile 01_structures/2myo_wt.pdb
set wt_all [atomselect top "all"]
set wt_q [vecsum [$wt_all get charge]]
puts "\n=========================================="
puts "VERIFICA TOPOLOGIA 2myo_wt (SOLO PROTEINA)"
puts "=========================================="
puts [format "Atomi totali ricostruiti (con H CHARMM36m): %d" [$wt_all num]]
puts [format "Carica netta della Myotrophin WT a pH fisiologico: %.2f e" $wt_q]
mol delete top

# 4. Solvatazione in box d'acqua esplicita TIP3P (padding = 12.0 A)
package require solvate
solvate 01_structures/2myo_wt.psf 01_structures/2myo_wt.pdb -t 12 -o 01_structures/2myo_wb

# 5. Neutralizzazione e aggiunta di 0.15 M NaCl con autoionize
package require autoionize
autoionize -psf 01_structures/2myo_wb.psf -pdb 01_structures/2myo_wb.pdb -sc 0.15 -o 01_structures/2myo_wb_ionized

# 6. Misura delle dimensioni esatte della cella periodica (PBC) per NAMD
mol delete all
mol new 01_structures/2myo_wb_ionized.psf
mol addfile 01_structures/2myo_wb_ionized.pdb

set all_sys [atomselect top "all"]
set wat_sel [atomselect top "water"]
set sod_sel [atomselect top "resname SOD"]
set cla_sel [atomselect top "resname CLA"]
set tot_q   [vecsum [$all_sys get charge]]

set minmax [measure minmax $all_sys]
set min_vec [lindex $minmax 0]
set max_vec [lindex $minmax 1]
set box_vec [vecsub $max_vec $min_vec]
set center  [measure center $all_sys]

set lx [lindex $box_vec 0]
set ly [lindex $box_vec 1]
set lz [lindex $box_vec 2]
set cx [lindex $center 0]
set cy [lindex $center 1]
set cz [lindex $center 2]

puts "\n=========================================="
puts "RIEPILOGO SISTEMA SOLVATATO E IONIZZATO (2myo_wb_ionized)"
puts "=========================================="
puts [format "Numero totale di atomi nel sistema:  %d" [$all_sys num]]
puts [format "Numero di molecole d'acqua (TIP3P):  %d" [expr {[$wat_sel num] / 3}]]
puts [format "Numero di ioni Na+ (SOD):            %d" [$sod_sel num]]
puts [format "Numero di ioni Cl- (CLA):            %d" [$cla_sel num]]
puts [format "Carica elettrica totale finale:      %.4f e" $tot_q]
puts "------------------------------------------"
puts "PARAMETRI PBC DA INSERIRE IN NAMD (.conf):"
puts [format "cellBasisVector1    %.3f    0.000    0.000" $lx]
puts [format "cellBasisVector2     0.000   %.3f    0.000" $ly]
puts [format "cellBasisVector3     0.000    0.000   %.3f" $lz]
puts [format "cellOrigin          %.3f   %.3f   %.3f" $cx $cy $cz]
puts "==========================================\n"

# Salva i parametri PBC anche su file per riferimento
set fp [open "01_structures/pbc_box_params.txt" w]
puts $fp "# Parametri PBC calcolati su 01_structures/2myo_wb_ionized"
puts $fp [format "cellBasisVector1    %.3f    0.000    0.000" $lx]
puts $fp [format "cellBasisVector2     0.000   %.3f    0.000" $ly]
puts $fp [format "cellBasisVector3     0.000    0.000   %.3f" $lz]
puts $fp [format "cellOrigin          %.3f   %.3f   %.3f" $cx $cy $cz]
close $fp

exit
