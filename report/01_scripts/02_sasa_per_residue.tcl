
# ==================================================
# 2MYO - Per-residue SASA Analysis
# ==================================================

# Get the currently loaded molecule
set mid [molinfo top]

if {$mid == -1} {
    error "No molecule loaded. Please load 2MYO.pdb first."
}

# Select protein heavy atoms and alpha carbons
set prot [atomselect $mid "protein and chain A and noh"]
set ca [atomselect $mid "protein and chain A and name CA"]

puts "Protein atoms: [$prot num]"
puts "Residues: [$ca num]"

# Store results before creating the CSV
set results {}
set total_sasa 0.0

# Calculate per-residue SASA
foreach residue [$ca get {resid resname}] {

    lassign $residue resid resname

    set sel [atomselect $mid \
        "protein and chain A and resid $resid and noh"]

    # Calculate SASA using a 1.4 A probe
    set status [catch {
        measure sasa 1.4 $prot -restrict $sel -samples 500
    } sasa]

    $sel delete

    if {$status != 0} {
        $prot delete
        $ca delete
        error "SASA calculation failed for residue $resid: $sasa"
    }

    # Print results in the console
    puts [format "%3d %-3s %.2f A^2" \
        $resid $resname $sasa]

    # Store results with three decimal places
    lappend results [format "%d,%s,%.3f" \
        $resid $resname $sasa]

    set total_sasa [expr {$total_sasa + $sasa}]
}

puts "Calculated residues: [llength $results]"
puts [format "Sum of per-residue SASA: %.3f A^2" $total_sasa]

# ==================================================
# CSV Export
# ==================================================

# Find the directory containing this TCL script
set script_dir [file dirname [file normalize [info script]]]

# Save results in notes/01_other/
set output_dir [file join $script_dir .. 01_other]

# Create the directory if it does not exist
file mkdir $output_dir

set filename [file normalize \
    [file join $output_dir "sasa_2myo_per_residue.csv"]]

# Open CSV file
set output [open $filename w]

# Write header
puts $output "Residue,ResName,SASA_A2"

# Write one row per residue
foreach row $results {
    puts $output $row
}

# Close file
close $output

# Confirm successful export
puts "========================================"
puts "CSV successfully saved."
puts "File: $filename"
puts "Residues exported: [llength $results]"
puts "File size: [file size $filename] bytes"
puts "========================================"

# Clean up atom selections
$prot delete
$ca delete