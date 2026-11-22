
# WT preparation verification — 2MYO
# Run from the BML repository root.

set outfile [open "report/02_utilities/wt_preparation_checks.txt" w]

proc report {line} {
    global outfile
    puts $line
    puts $outfile $line
}

# Load the prepared WT system
set mid [mol new "01_structures/2myo_wb_ionized.psf" type psf waitfor all]
mol addfile "01_structures/2myo_wb_ionized.pdb" type pdb waitfor all molid $mid

set all [atomselect $mid "all"]
set protein [atomselect $mid "segname A"]
set heavy [atomselect $mid "segname A and noh"]
set ca [atomselect $mid "segname A and name CA"]
set waters [atomselect $mid "water"]
set sodium [atomselect $mid "resname SOD"]
set chloride [atomselect $mid "resname CLA"]

# 1. Composition and charges
report "=== SYSTEM COMPOSITION ==="
report "Total atoms: [$all num]"
report "Protein atoms: [$protein num]"
report "Protein heavy atoms: [$heavy num]"
report "CA atoms: [$ca num]"
report "Water molecules: [expr {[$waters num] / 3}]"
report "Na+: [$sodium num]"
report "Cl-: [$chloride num]"
report [format "Protein charge: %.4f e" [vecsum [$protein get charge]]]
report [format "Total charge: %.4f e" [vecsum [$all get charge]]]

# 2. Histidine protonation
report "\n=== HISTIDINE STATES ==="

foreach resid {40 67 68 80} {
    set h [atomselect $mid "segname A and resid $resid"]
    set hd1 [atomselect $mid "segname A and resid $resid and name HD1"]
    set he2 [atomselect $mid "segname A and resid $resid and name HE2"]

    report "Residue $resid: [lsort -unique [$h get resname]], HD1=[$hd1 num], HE2=[$he2 num]"

    $h delete
    $hd1 delete
    $he2 delete
}

# 3. Structural distances relevant to histidines
report "\n=== HISTIDINE CONTACT DISTANCES ==="

foreach contact {{40 ND1 74 OG} {67 ND1 37 NZ} {67 NE2 35 O}} {
    lassign $contact rA aA rB aB

    set sa [atomselect $mid "segname A and resid $rA and name $aA"]
    set sb [atomselect $mid "segname A and resid $rB and name $aB"]

    if {[$sa num] == 1 && [$sb num] == 1} {
        set a [lindex [$sa get {x y z}] 0]
        set b [lindex [$sb get {x y z}] 0]
        set d [veclength [vecsub $a $b]]
        report [format "%s%d - %s%d: %.3f A" $aA $rA $aB $rB $d]
    }

    $sa delete
    $sb delete
}

# 4. Sequence and coordinate preservation
report "\n=== STRUCTURAL INTEGRITY ==="

set rawid [mol new "01_structures/raw/2MYO.pdb" type pdb waitfor all]
set rawca [atomselect $rawid "chain A and name CA"]

set originalResidues [$rawca get resid]
set preparedResidues [$ca get resid]

if {$originalResidues eq $preparedResidues} {
    set maxdist 0.0

    foreach a [$rawca get {x y z}] b [$ca get {x y z}] {
        set d [veclength [vecsub $a $b]]
        if {$d > $maxdist} {
            set maxdist $d
        }
    }

    report "Matching CA residues: [llength $originalResidues]"
    report [format "Maximum CA displacement: %.4f A" $maxdist]
} else {
    report "WARNING: residue numbering or order does not match"
}

# 5. Solvation box
report "\n=== SOLVATION BOX ==="

set mm [measure minmax $all]
set mn [lindex $mm 0]
set mx [lindex $mm 1]

report "Box dimensions (A): [vecsub $mx $mn]"
report "Atomic-coordinate center (A): [measure center $all]"
report "Bounding-box midpoint (A): [vecscale 0.5 [vecadd $mn $mx]]"

set pm [measure minmax $heavy]

report "Minimum-side protein margins (A): [vecsub [lindex $pm 0] $mn]"
report "Maximum-side protein margins (A): [vecsub $mx [lindex $pm 1]]"

# 6. Close-contact screening
report "\n=== CLOSE CONTACTS ==="

set contactPairs [measure contacts 2.25 $heavy]
set allres [$all get resid]
set allnames [$all get name]
set allcoords [$all get {x y z}]

set count 0

foreach i [lindex $contactPairs 0] j [lindex $contactPairs 1] {
    set ri [lindex $allres $i]
    set rj [lindex $allres $j]

    if {[expr {abs($ri - $rj)}] > 1} {
        incr count
        set d [veclength [vecsub [lindex $allcoords $i] [lindex $allcoords $j]]]

        report [format "%s%d - %s%d: %.3f A" \
            [lindex $allnames $i] $ri \
            [lindex $allnames $j] $rj $d]
    }
}

report "Close heavy-atom pairs (<2.25 A): $count"

set waterO [atomselect $mid "water and name OH2"]
set ions [atomselect $mid "resname SOD CLA"]

report "Protein-water O pairs (<2.0 A): [llength [lindex [measure contacts 2.0 $heavy $waterO] 0]]"
report "Protein-ion pairs (<2.0 A): [llength [lindex [measure contacts 2.0 $heavy $ions] 0]]"

close $outfile

puts "\nResults saved in 02_utilities/wt_preparation_checks.txt"