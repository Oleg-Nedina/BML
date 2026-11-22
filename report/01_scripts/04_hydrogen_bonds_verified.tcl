
# ==================================================
# 2MYO - N-H...O Hydrogen Bond Analysis
# ==================================================

set mid [molinfo top]

if {$mid < 0} {
    error "Load the original 2MYO.pdb first."
}

# Separate nitrogen donors and oxygen acceptors
set donors [atomselect $mid \
    "protein and chain A and element N"]

set acceptors [atomselect $mid \
    "protein and chain A and element O"]

# Atom metadata indexed by VMD atom index
set all [atomselect $mid "all"]
set info [dict create]

foreach atom [$all get {index resid resname name}] {
    lassign $atom idx resid resname name
    dict set info $idx [list $resid $resname $name]
}

# Identify geometric hydrogen-bond candidates
set hb [measure hbonds 3.5 30 $donors $acceptors]

set dlist [lindex $hb 0]
set alist [lindex $hb 1]
set hlist [lindex $hb 2]

# Output directory
set script_dir [file dirname [file normalize [info script]]]
set output_dir [file join $script_dir .. 01_utilities]
file mkdir $output_dir

set outfile [file normalize \
    [file join $output_dir "2myo_hbonds_verified.csv"]]

set f [open $outfile w]
puts $f "DonorResidue,DonorAtom,AcceptorResidue,AcceptorAtom,DA_Distance_A,DHA_Angle_deg"

set seen [dict create]
set count 0

foreach d $dlist a $alist h $hlist {

    if {$h < 0} {
        continue
    }

    lassign [dict get $info $d] dr dn da
    lassign [dict get $info $a] ar an aa

    # Exclude intra-residue interactions
    if {$dr == $ar} {
        continue
    }

    # Enforce chemically plausible atom types
    if {![string match "N*" $da]} {
        continue
    }

    if {![string match "O*" $aa]} {
        continue
    }

    # Verify donor-hydrogen geometry
    set dh [measure bond [list $d $h] molid $mid]

    if {$dh < 0.65 || $dh > 1.25} {
        continue
    }

    # Verify donor-acceptor distance and D-H-A angle
    set dist [measure bond [list $d $a] molid $mid]
    set angle [measure angle [list $d $h $a] molid $mid]

    if {$dist > 3.5 || $angle < 150} {
        continue
    }

    # Avoid duplicate interactions
    set key "$d-$a-$h"

    if {[dict exists $seen $key]} {
        continue
    }

    dict set seen $key 1

    puts $f [format "%s%d,%s,%s%d,%s,%.3f,%.2f" \
        $dn $dr $da $an $ar $aa $dist $angle]

    incr count
}

close $f

$donors delete
$acceptors delete
$all delete

puts "===================================="
puts "N-H...O hydrogen bonds: $count"
puts "CSV saved: $outfile"
puts "===================================="