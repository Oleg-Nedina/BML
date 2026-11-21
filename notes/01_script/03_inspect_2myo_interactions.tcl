
# ============================================================
# 2MYO - Static hydrogen bonds and candidate salt bridges
# ============================================================
# Run in VMD Tk Console AFTER loading the original 2MYO.pdb.
# This script analyzes ONLY the current frame of the top molecule.
# No molecular-dynamics trajectory is required.
#
# Outputs:
# - 2myo_hbonds_static.csv
# - 2myo_saltbridges_static.csv
# ============================================================

# STEP 1 - Select the protein

set mid [molinfo top]

if {$mid < 0} {
    error "Load 2MYO.pdb in VMD before running this script."
}

set protein [atomselect $mid "protein and chain A"]

if {[$protein num] == 0} {
    error "No protein atoms found in chain A of the top molecule."
}

# Fetch atom metadata.
# The "all" selection is ordered by atom index.

set allAtoms [atomselect $mid "all"]
set info [$allAtoms get {resid resname name x y z}]


# ============================================================
# STEP 2 - HYDROGEN BOND ANALYSIS
# ============================================================
#
# Scientific question:
# Which residues form geometrically compatible hydrogen bonds?
#
# Criteria:
# - Donor-acceptor distance < 3.5 Angstrom
# - Angular deviation from linear D-H-A < 30 degrees
#
# These are geometric candidates, not interaction energies.
# ============================================================

set hbonds [measure hbonds 3.5 30 $protein]

set donors    [lindex $hbonds 0]
set acceptors [lindex $hbonds 1]
set hydrogens [lindex $hbonds 2]

# Create output CSV

set hFile "2myo_hbonds_static.csv"
set hout [open $hFile w]

puts $hout "DonorResidue,DonorAtom,AcceptorResidue,AcceptorAtom,HydrogenAtom,DA_Distance_A"

set nHB 0
set seenHB [dict create]

foreach d $donors a $acceptors h $hydrogens {

    # Retrieve donor information
    lassign [lindex $info $d] dr dn da dx dy dz

    # Retrieve acceptor information
    lassign [lindex $info $a] ar an aa ax ay az

    # Exclude interactions within the same residue
    if {$dr eq $ar} {
        continue
    }

    # Avoid duplicate hydrogen bonds
    set hbkey "$d-$a-$h"

    if {[dict exists $seenHB $hbkey]} {
        continue
    }

    dict set seenHB $hbkey 1

    # Retrieve hydrogen atom name
    set hname "unknown"

    if {$h >= 0 && $h < [llength $info]} {
        set hname [lindex [lindex $info $h] 2]
    }

    # Calculate donor-acceptor distance
    set dist [veclength [vecsub \
        [list $dx $dy $dz] \
        [list $ax $ay $az]]]

    # Write result to CSV
    puts $hout [format "%s%s,%s,%s%s,%s,%s,%.3f" \
        $dn $dr $da $an $ar $aa $hname $dist]

    incr nHB
}

close $hout


# ============================================================
# STEP 3 - SALT BRIDGE ANALYSIS
# ============================================================
#
# Scientific question:
# Which oppositely charged residues form close contacts?
#
# Acidic residues:
# - ASP: OD1, OD2
# - GLU: OE1, OE2
#
# Basic residues:
# - LYS: NZ
# - ARG: NE, NH1, NH2
#
# Distance cutoff: 4.0 Angstrom
#
# Histidines are excluded because their protonation
# states require further investigation.
# ============================================================

# Select acidic side-chain oxygen atoms

set acidic [atomselect $mid \
    "protein and chain A and ((resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2))"]

# Select basic side-chain nitrogen atoms

set basic [atomselect $mid \
    "protein and chain A and ((resname LYS and name NZ) or (resname ARG and name NE NH1 NH2))"]

# Identify contacts within 4 Angstrom

set contacts [measure contacts 4.0 $acidic $basic]

set negIndices [lindex $contacts 0]
set posIndices [lindex $contacts 1]

# Store the shortest distance for each residue pair

set best [dict create]

foreach neg $negIndices pos $posIndices {

    # Retrieve acidic atom information
    lassign [lindex $info $neg] nr nn na nx ny nz

    # Retrieve basic atom information
    lassign [lindex $info $pos] pr pn pa px py pz

    # Ignore contacts within the same residue
    if {$nr eq $pr} {
        continue
    }

    # Calculate the distance
    set dist [veclength [vecsub \
        [list $nx $ny $nz] \
        [list $px $py $pz]]]

    # Create a unique identifier for each residue pair
    set key "${nn}${nr}-${pn}${pr}"

    # Keep only the shortest atom-to-atom distance
    if {![dict exists $best $key] || \
        $dist < [lindex [dict get $best $key] 0]} {

        dict set best $key \
            [list $dist $nn $nr $na $pn $pr $pa]
    }
}

# Create salt bridge CSV

set sFile "2myo_saltbridges_static.csv"
set sout [open $sFile w]

puts $sout "AcidicResidue,AcidicAtom,BasicResidue,BasicAtom,Minimum_ON_Distance_A"

foreach key [lsort [dict keys $best]] {

    lassign [dict get $best $key] dist nn nr na pn pr pa

    puts $sout [format "%s%s,%s,%s%s,%s,%.3f" \
        $nn $nr $na $pn $pr $pa $dist]
}

close $sout


# ============================================================
# STEP 4 - PRINT RESULTS
# ============================================================

puts "\n============================================"
puts "2MYO - STATIC INTERACTION ANALYSIS"
puts "============================================"

puts "\nHYDROGEN BONDS"
puts "Inter-residue H-bond candidates: $nHB"

puts "\nSALT BRIDGES"
puts "Candidate salt-bridge residue pairs: [dict size $best]"

foreach key [lsort [dict keys $best]] {

    lassign [dict get $best $key] dist nn nr na pn pr pa

    puts [format "  %-15s %s-%s   %.2f A" \
        $key $na $pa $dist]
}

puts "\nOUTPUT FILES"
puts "Hydrogen bonds: $hFile"
puts "Salt bridges:   $sFile"
puts "Saved in: [pwd]"

puts "\nIMPORTANT NOTES:"
puts "- Results are based on one NMR structure."
puts "- Hydrogen bonds are geometric candidates."
puts "- Salt bridges are close ionic contacts."
puts "- Interaction persistence is not evaluated."
puts "- Stability contributions cannot be inferred directly."

puts "\nAnalysis completed."
puts "============================================"


# ============================================================
# STEP 5 - CLEAN UP
# ============================================================

$protein delete
$acidic delete
$basic delete
$allAtoms delete
