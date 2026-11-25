# ==============================================================================
# viz_05_md_comparison_400K.tcl
# Visualizzazione Dinamica 3D in VMD della traiettoria MD a 400 K (2.5 ns, 1252 frame)
# del Triplo Mutante vincitore (K37L_H68G_T103L) confrontato con il Wild-Type (WT).
#
# Funzionalità automatiche:
#   - Carica il Triplo Mutante K37L_H68G_T103L a 400 K (e opzionalmente il WT a 400 K)
#   - Nasconde l'acqua per una vista cristallina della proteina
#   - Allinea tutti i 1252 frame della traiettoria sul backbone del Frame 0
#   - Colora i 4 Ankyrin Repeats e mostra in Licorice spesso i 3 residui mutati
#     (LEU 37 in ANK2, GLY 68 in ANK3, LEU 103 in ANK4) + ASP 60 (la cerniera stabilizzata!)
#
# Esecuzione dalla cartella BML/:
#   vmd -e scripts/viz_05_md_comparison_400K.tcl
# (Oppure per caricare un altro mutante singolo/doppio:
#   vmd -e scripts/viz_05_md_comparison_400K.tcl -args K37L )
# ==============================================================================

color Display Background white
axes location Off
display projection Orthographic
display depthcue off
display rendermode GLSL

set target_mut "K37L_H68G_T103L"
set seed_num   1
if {$argc >= 1} {
    set target_mut [lindex $argv 0]
}
if {$argc >= 2} {
    set seed_num [lindex $argv 1]
}

if {$target_mut eq "WT"} {
    set mut_psf "02_wt_baseline/wt_wb_ionized.psf"
    set mut_dcd "04_md_simulations/WT/wt_400K_seed${seed_num}.dcd"
    if {![file exists $mut_dcd] && $seed_num == 1} {
        set mut_dcd "02_wt_baseline/wt_400K.dcd"
    }
} else {
    set mut_psf "04_md_simulations/${target_mut}/${target_mut}_wb_ionized.psf"
    set mut_dcd "04_md_simulations/${target_mut}/${target_mut}_400K_seed${seed_num}.dcd"
}

if {![file exists $mut_dcd]} {
    puts "ERRORE: Traiettoria $mut_dcd non trovata!"
    return
}

# 1. Carica la traiettoria a 400 K del Mutante selezionato
mol new $mut_psf
mol addfile $mut_dcd step 2 waitfor all
mol rename top "${target_mut}_400K"
mol delrep 0 top

# Allinea tutti i frame sul backbone della proteina al Frame 0
set nf [molinfo top get numframes]
set all_sel [atomselect top "all"]
set bb_ref  [atomselect top "protein and backbone" frame 0]
set bb_cur  [atomselect top "protein and backbone"]
for {set f 0} {$f < $nf} {incr f} {
    $all_sel frame $f
    $bb_cur frame $f
    $all_sel move [measure fit $bb_cur $bb_ref]
}
$all_sel delete; $bb_ref delete; $bb_cur delete

# Rappresentazione 1: le 4 regioni strutturali Yang et al. (1998)
# (N-term=Blu, Repeat 1=Ciano, Repeat 2=Arancione, C-term=Verde)
foreach {r_range cid} {"1 to 27" 0 "28 to 60" 10 "61 to 93" 3 "94 to 118" 7} {
    mol representation NewCartoon 0.30 20.0 4.1 0
    mol color ColorID $cid
    mol selection "protein and resid $r_range"
    mol material Opaque
    mol addrep top
}

# Rappresentazione 2: I 3 Siti Mutati (LEU 37, GLY 68, LEU 103) in Licorice spesso
mol representation Licorice 0.32 20.0 20.0
mol color Name
mol selection "protein and resid 37 68 103 and noh"
mol material Opaque
mol addrep top

# Rappresentazione 3: Il core idrofobico di ANK2 salvato da K37L (ILE 50, LEU 51, PHE 53) e la cerniera ASP 60
mol representation VDW 0.45 16.0
mol color ColorID 4
mol selection "protein and resid 37 50 51 53 60 68 102 103 and name CA"
mol material Opaque
mol addrep top

# Rappresentazione 4: Sito Funzionale 3AAA (24 residui che legano Actin-Capping Protein) in superficie trasparente
mol representation QuickSurf 0.9 0.5 1.0 1.0
mol color ColorID 10
mol selection "protein and resid 4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100"
mol material Glass1
mol addrep top

display resetview
animate goto 0

puts "\n=========================================================================="
puts " TRAIETTORIA MD A 400 K CARICATA E ALLINEATA: ${target_mut} Seed ${seed_num} ($nf frame)"
puts "   - Blu      : N-terminal (1-27, H1-H2)       | Ciano : Repeat 1 (28-60, H3-H4, con K37L)"
puts "   - Arancione: Repeat 2 (61-93, H5-H6, H68G)  | Verde : C-terminal (94-118, H7, con T103L)"
puts "   - Licorice spesso: Catene laterali nei siti 37, 68, 103"
puts "   - Sfere Giallo/Oro: C-alpha dei siti chiave (37, 50, 51, 53, 60, 68, 102, 103)"
puts "   - Superficie trasparente Ciano: Sito Funzionale 3AAA (24 residui)"
puts " Uso opzionale: vmd -e scripts/viz_05_md_comparison_400K.tcl -args <WT|K37L|T103L|K37L_T103L|K37L_H68G_T103L> <seed>"
puts " Premi il tasto PLAY in basso a destra nella finestra VMD Main per vedere il film!"
puts "==========================================================================\n"
