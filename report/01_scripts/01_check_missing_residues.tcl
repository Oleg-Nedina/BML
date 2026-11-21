
# ============================================================
# 2MYO - Protein Structural Completeness Analysis
# ============================================================
# Objective:
# 1. Compare SEQRES with ATOM records.
# 2. Identify missing amino acid residues.
# 3. Identify missing backbone and side-chain heavy atoms.
# 4. Detect unexpected or duplicate heavy atoms.
# 5. Check for alternative atom locations.
#
# Hydrogen completeness is NOT evaluated because it depends
# on protonation states and molecular preparation.
# ============================================================

# STEP 1 - Select the original PDB file

set pdbFile [tk_getOpenFile \
    -title "Select the original 2MYO.pdb file" \
    -filetypes {{{PDB files} {.pdb}} {{All files} {*}}}]

if {$pdbFile eq "" || ![file exists $pdbFile]} {
    error "No valid PDB file selected."
}

# STEP 2 - Define expected heavy atoms for each amino acid

set expected [dict create \
    ALA {N CA C O CB} \
    ARG {N CA C O CB CG CD NE CZ NH1 NH2} \
    ASN {N CA C O CB CG OD1 ND2} \
    ASP {N CA C O CB CG OD1 OD2} \
    CYS {N CA C O CB SG} \
    GLN {N CA C O CB CG CD OE1 NE2} \
    GLU {N CA C O CB CG CD OE1 OE2} \
    GLY {N CA C O} \
    HIS {N CA C O CB CG ND1 CD2 CE1 NE2} \
    ILE {N CA C O CB CG1 CG2 CD1} \
    LEU {N CA C O CB CG CD1 CD2} \
    LYS {N CA C O CB CG CD CE NZ} \
    MET {N CA C O CB CG SD CE} \
    PHE {N CA C O CB CG CD1 CD2 CE1 CE2 CZ} \
    PRO {N CA C O CB CG CD} \
    SER {N CA C O CB OG} \
    THR {N CA C O CB OG1 CG2} \
    TRP {N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2} \
    TYR {N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH} \
    VAL {N CA C O CB CG1 CG2}]

# STEP 3 - Read SEQRES and ATOM records

set seqres {}
set observed [dict create]
set residueTypes [dict create]

set nHeavy 0
set nHydrogen 0
set nCA 0

set extraChains {}
set insertions {}
set altlocs {}

set f [open $pdbFile r]

while {[gets $f line] >= 0} {

    set record [string trim [string range $line 0 5]]

    # Read expected amino acid sequence from SEQRES

    if {$record eq "SEQRES" &&
        [string trim [string range $line 11 11]] eq "A"} {

        foreach aa [regexp -all -inline {\S+} \
            [string range $line 19 end]] {
            lappend seqres $aa
        }
    }

    # Analyze atomic coordinates

    if {$record ne "ATOM"} {
        continue
    }

    set chain [string trim [string range $line 21 21]]

    if {$chain ne "A"} {
        lappend extraChains $chain
        continue
    }

    set resid [string trim [string range $line 22 25]]
    set resname [string trim [string range $line 17 19]]
    set atom [string trim [string range $line 12 15]]
    set element [string toupper \
        [string trim [string range $line 76 77]]]

    set icode [string trim [string range $line 26 26]]
    set altloc [string trim [string range $line 16 16]]

    if {$icode ne ""} {
        lappend insertions "$resid$icode"
    }

    if {$altloc ne ""} {
        lappend altlocs [format "%s%s:%s(%s)" \
            $resname $resid $atom $altloc]
    }

    # Count hydrogen atoms separately

    if {$element eq "H" || $element eq "D"} {
        incr nHydrogen
        continue
    }

    incr nHeavy

    if {$atom eq "CA"} {
        incr nCA
    }

    dict lappend observed $resid $atom
    dict set residueTypes $resid $resname
}

close $f

# STEP 4 - Compare expected and observed residues

set missingResidues {}
set mismatchedResidues {}
set missingAtoms {}
set unexpectedAtoms {}
set duplicateAtoms {}
set unsupportedResidues {}

set nExpected [llength $seqres]

# Specific to 2MYO: continuous residue numbering 1-118

for {set i 1} {$i <= $nExpected} {incr i} {

    set expectedAA [lindex $seqres [expr {$i - 1}]]

    # Check whether the residue exists

    if {![dict exists $observed $i]} {
        lappend missingResidues "$expectedAA$i"
        continue
    }

    set actualAA [dict get $residueTypes $i]

    # Check amino acid identity

    if {$actualAA ne $expectedAA} {
        lappend mismatchedResidues \
            "$i (SEQRES=$expectedAA, ATOM=$actualAA)"
    }

    if {![dict exists $expected $expectedAA]} {
        lappend unsupportedResidues "$expectedAA$i"
        continue
    }

    # Expected atoms for this amino acid

    set required [dict get $expected $expectedAA]

    # Check terminal oxygen in the last residue

    if {$i == $nExpected} {
        lappend required OXT
    }

    set found [dict get $observed $i]

    # Check for missing heavy atoms

    foreach atom $required {
        if {[lsearch -exact $found $atom] == -1} {
            lappend missingAtoms "$expectedAA$i:$atom"
        }
    }

    # Check for unexpected and duplicate atoms

    foreach atom [lsort -unique $found] {

        if {[lsearch -exact $required $atom] == -1} {
            lappend unexpectedAtoms "$actualAA$i:$atom"
        }

        if {[llength [lsearch -all -exact $found $atom]] > 1} {
            lappend duplicateAtoms "$actualAA$i:$atom"
        }
    }
}

# STEP 5 - Check for residues outside expected numbering

set extraResidues {}

foreach i [dict keys $observed] {
    if {![string is integer -strict $i] ||
        $i < 1 || $i > $nExpected} {
        lappend extraResidues $i
    }
}

# STEP 6 - Print results

puts "\n========== 2MYO STRUCTURAL COMPLETENESS CHECK =========="
puts "PDB file: $pdbFile"
puts "Expected residues (SEQRES, chain A): $nExpected"
puts "Residues with heavy-atom coordinates: [dict size $observed]"
puts "Number of CA atoms: $nCA"
puts "Number of heavy atoms: $nHeavy"
puts "Number of H/D atoms: $nHydrogen"

foreach {title values} [list \
    "Missing residues" $missingResidues \
    "Residue identity mismatches" $mismatchedResidues \
    "Missing heavy atoms" $missingAtoms \
    "Unexpected heavy atoms" $unexpectedAtoms \
    "Duplicate heavy atoms" $duplicateAtoms \
    "Residues outside the expected range" $extraResidues \
    "Non-standard residues (manual check)" $unsupportedResidues \
    "Insertion codes" $insertions \
    "Alternate atom locations" $altlocs \
    "Other ATOM chains" $extraChains] {

    if {[llength $values] == 0} {
        puts "$title: NONE"
    } else {
        puts "$title: [join [lsort -unique $values] {, }]"
    }
}

# STEP 7 - Final assessment

if {$nExpected > 0 &&
    $nCA == $nExpected &&
    [llength $missingResidues] == 0 &&
    [llength $mismatchedResidues] == 0 &&
    [llength $missingAtoms] == 0 &&
    [llength $unexpectedAtoms] == 0 &&
    [llength $duplicateAtoms] == 0 &&
    [llength $extraResidues] == 0 &&
    [llength $unsupportedResidues] == 0 &&
    [llength $insertions] == 0 &&
    [llength $altlocs] == 0} {

    puts "RESULT: COMPLETE standard residues and heavy atoms."

} else {

    puts "RESULT: CHECK the reported issues before preparing MD."
}

puts "NOTE: Hydrogen protonation, bond geometry, and force-field compatibility are NOT checked."
puts "=======================================================\n"
