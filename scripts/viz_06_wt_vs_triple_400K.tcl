# ==============================================================================
# viz_06_wt_vs_triple_400K.tcl
# Sovrapposizione Dinamica 3D in VMD: Wild-Type (WT) vs Triplo Mutante
# (K37L_H68G_T103L) a 400 K (2.5 ns), allineati sul backbone del Frame 0 del WT!
#
# Ideale per fare screenshot comparativi ad alto impatto per le slide finali e
# per osservare in tempo reale:
#   - L'espansione/unfolding del WT (Rosso semitrasparente) vs la compattezza
#     del Triplo Mutante (Viola/Indaco solido)
#   - Il comportamento delle ultime eliche H6 (82-89) e H7 (112-117, C-terminale)
#   - I 3 siti mutati (LEU 37, GLY 68, LEU 103) e la cerniera ASP 60
#   - L'integrità della Functional Mask 3AAA (24 residui, toggleabile)
#
# Esecuzione dalla cartella BML/:
#   vmd -e scripts/viz_06_wt_vs_triple_400K.tcl
# (Oppure scegliendo un seed specifico 1, 2 o 3:
#   vmd -e scripts/viz_06_wt_vs_triple_400K.tcl -args 2 )
# ==============================================================================

color Display Background white
axes location Off
display projection Orthographic
display depthcue off
display rendermode GLSL
display antialias on
display ambientocclusion on

set seed_num 1
if {$argc >= 1} {
    set seed_num [lindex $argv 0]
}

# Percorsi WT (compatibile sia con 04_md_simulations/WT/ sia con 02_wt_baseline/ per Seed 1)
set wt_psf "02_wt_baseline/wt_wb_ionized.psf"
set wt_dcd "04_md_simulations/WT/wt_400K_seed${seed_num}.dcd"
if {![file exists $wt_dcd] && $seed_num == 1} {
    set wt_dcd "02_wt_baseline/wt_400K.dcd"
}

set tr_psf "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_wb_ionized.psf"
set tr_dcd "04_md_simulations/K37L_H68G_T103L/K37L_H68G_T103L_400K_seed${seed_num}.dcd"

if {![file exists $wt_dcd]} {
    puts "ERRORE: Traiettoria WT non trovata ($wt_dcd)!"
    return
}
if {![file exists $tr_dcd]} {
    puts "ERRORE: Traiettoria Triplo Mutante non trovata ($tr_dcd)!"
    return
}

# 1. Carica Wild-Type (Mol 0)
set mol_wt [mol new $wt_psf]
mol addfile $wt_dcd step 2 waitfor all molid $mol_wt
mol rename $mol_wt "WT_400K_seed${seed_num}"
mol delrep 0 $mol_wt

# 2. Carica Triplo Mutante (Mol 1)
set mol_tr [mol new $tr_psf]
mol addfile $tr_dcd step 2 waitfor all molid $mol_tr
mol rename $mol_tr "Triple_K37L_H68G_T103L_400K_seed${seed_num}"
mol delrep 0 $mol_tr

# 3. Allinea TUTTI i frame di entrambe le molecole sul backbone (CA) del Frame 0 del WT
set ref_ca [atomselect $mol_wt "protein and name CA" frame 0]

foreach mid [list $mol_wt $mol_tr] {
    set nf [molinfo $mid get numframes]
    set all_sel [atomselect $mid "all"]
    set cur_ca  [atomselect $mid "protein and name CA"]
    for {set f 0} {$f < $nf} {incr f} {
        $all_sel frame $f
        $cur_ca  frame $f
        $all_sel move [measure fit $cur_ca $ref_ca]
    }
    $all_sel delete
    $cur_ca delete
}
$ref_ca delete

# --- RAPPRESENTAZIONI WILD-TYPE (Mol 0: Rosso Trasparente) ---
# Rep 0: Backbone completo WT in Rosso semitrasparente (Ghost)
mol representation NewCartoon 0.28 20.0 4.1 0
mol color ColorID 1
mol selection "protein"
mol material Transparent
mol addrep $mol_wt

# Rep 1: Ultime eliche WT — H6 (82-89) e H7 (112-117) in Rosso opaco (per vedere lo sfilacciamento!)
mol representation NewCartoon 0.34 20.0 4.1 0
mol color ColorID 1
mol selection "protein and (resid 82 to 89 or resid 112 to 117)"
mol material AOChalky
mol addrep $mol_wt
mol showrep $mol_wt 1 off

# --- RAPPRESENTAZIONI TRIPLO MUTANTE (Mol 1: Viola Solido + Dettagli) ---
# Rep 0: Backbone completo Triplo Mutante in Viola/Indaco (ColorID 11)
mol representation NewCartoon 0.32 20.0 4.1 0
mol color ColorID 11
mol selection "protein"
mol material AOChalky
mol addrep $mol_tr

# Rep 1: Ultime eliche Triplo — H6 (82-89) e H7 (112-117) in Verde brillante (ColorID 7)
mol representation NewCartoon 0.36 20.0 4.1 0
mol color ColorID 7
mol selection "protein and (resid 82 to 89 or resid 112 to 117)"
mol material AOChalky
mol addrep $mol_tr
mol showrep $mol_tr 1 off

# Rep 2: I 3 siti mutati (LEU 37, GLY 68, LEU 103) in Licorice spesso (colori atomici)
mol representation Licorice 0.30 20.0 20.0
mol color Name
mol selection "protein and resid 37 68 103 and noh"
mol material AOShiny
mol addrep $mol_tr

# Rep 3: Cerniera stabilizzata ASP 60 in sfere VDW Oro/Giallo (ColorID 4)
mol representation VDW 0.45 16.0
mol color ColorID 4
mol selection "protein and resid 60 and name CA"
mol material Opaque
mol addrep $mol_tr

# Rep 4: Functional Mask 3AAA (24 residui) in superficie trasparente Ciano (spenta di default)
mol representation QuickSurf 0.9 0.5 1.0 1.0
mol color ColorID 10
mol selection "protein and resid 4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100"
mol material Glass1
mol addrep $mol_tr
mol showrep $mol_tr 4 off

display resetview
animate goto last

puts "\n=========================================================================="
puts " CONFRONTO 3D VMD CARICATO: WT (Rosso) vs TRIPLO MUTANTE (Viola) — Seed ${seed_num}"
puts "   Posizionato sull'ULTIMO FRAME (2.5 ns a 400 K) per impatto visivo immediato!"
puts ""
puts "   - Rosso trasparente : Wild-Type (2MYO) a 400 K (espanso/deformato)"
puts "   - Viola solido      : Triplo Mutante K37L+H68G+T103L (compatto!)"
puts "   - Licorice spesso   : Siti mutati LEU 37 (Rep1), GLY 68 (Rep2), LEU 103 (C-term)"
puts "   - Sfera Gialla      : Cerniera ASP 60 (RMSF ridotto da 2.53/5.00 a 1.50 A)"
puts ""
puts "   SUGGERIMENTI IN 'Graphics -> Representations':"
puts "   - Accendi Rep 1 su WT e Triplo per evidenziare le ultime eliche H6 (82-89) e H7 (112-117)"
puts "   - Accendi Rep 4 sul Triplo per mostrare la Functional Mask 3AAA intatta"
puts "   - Premi PLAY in VMD Main per guardare le due traiettorie in sincronia!"
puts "   - Salva screenshot per le slide in: 06_presentation/presentazione_prof_WT_vs_Triplo/"
puts "==========================================================================\n"
