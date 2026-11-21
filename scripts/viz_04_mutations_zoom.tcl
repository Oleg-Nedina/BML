# ==============================================================================
# viz_04_mutations_zoom.tcl
# Punto 3 della Scaletta del Prof (Laboratory_Project_2026.pdf, Slide 3 & 19):
# Confronto strutturale 3D in VMD tra il Matched Wild-Type e il Triplo/Doppio
# Mutante generato da FoldX BuildModel (K37L in ANK2, H68G in ANK3, T103L in ANK4),
# evidenziando come le sostituzioni migliorano le interazioni locali senza
# toccare i 24 residui della Functional Mask (3AAA).
#
# Esecuzione dalla cartella principale BML/:
#   vmd -e scripts/viz_04_mutations_zoom.tcl
# ==============================================================================

color Display Background white
axes location Off
display projection Orthographic
display depthcue off
display rendermode GLSL

set wt_pdb  "03_mutations/final_models_pdb/2myo_matched_WT_for_K37L_H68G_T103L.pdb"
set mut_pdb "03_mutations/final_models_pdb/2myo_mut_K37L_H68G_T103L.pdb"

# 1. Carica il Matched Wild-Type (Mol 0 - Semitrasparente Grigio/Rosso sui siti WT)
mol new $wt_pdb
mol rename top "WT_Matched_Reference"
mol delrep 0 top

# Backbone WT trasparente
mol representation NewCartoon 0.28 20.0 4.1 0
mol color ColorID 2
mol selection "protein"
mol material Transparent
mol addrep top

# Catene laterali WT nei 3 siti mutati (K37, H68, T103) in Rosso/Arancio (Stato Originale Instabile)
mol representation Licorice 0.22 16.0 16.0
mol color ColorID 1
mol selection "resid 37 68 103 and noh"
mol material Opaque
mol addrep top

# 2. Carica il Modello Mutato BuildModel (Mol 1 - Colorato per Ankyrin Repeat + Siti Mutati in Verde/Blu)
mol new $mut_pdb
mol rename top "Mutant_K37L_H68G_T103L"
mol delrep 0 top

# ANK1 (1-34, Blu), ANK2 (35-67, Ciano), ANK3 (68-100, Arancio), ANK4 (101-118, Verde)
foreach {r_range cid} {"1 to 34" 0 "35 to 67" 10 "68 to 100" 3 "101 to 118" 7} {
    mol representation NewCartoon 0.30 20.0 4.1 0
    mol color ColorID $cid
    mol selection "resid $r_range"
    mol material Opaque
    mol addrep top
}

# Catene laterali MUTATE (LEU 37, GLY 68, LEU 103) in Verde/Name (Licorice spesso)
mol representation Licorice 0.30 20.0 20.0
mol color Name
mol selection "resid 37 68 103 and noh"
mol material Opaque
mol addrep top

# Intorno locale entro 4.5 A dai siti mutati (per mostrare il nuovo packing idrofobico e l'assenza di clash)
mol representation Licorice 0.14 14.0 14.0
mol color Name
mol selection "(noh and not resid 37 68 103) and within 4.5 of (resid 37 68 103)"
mol material Transparent
mol addrep top

# Ponte H inter-repeat protetto: HIS 40 (HSD) -> SER 74 (OG)
mol representation Licorice 0.18 16.0 16.0
mol color ColorID 11
mol selection "resid 40 74 and noh"
mol material Opaque
mol addrep top

display resetview
puts "\n=========================================================================="
puts " SCENA 3D VMD CARICATA: Confronto Matched WT (Rosso) vs Mutanti (Name)"
puts "   - Residuo 37  (ANK2): LYS 37 (WT Rosso) -> LEU 37 (Mutato)"
puts "   - Residuo 68  (ANK3): HIS 68 (WT Rosso) -> GLY 68 (Mutato)"
puts "   - Residuo 103 (ANK4): THR 103 (WT Rosso) -> LEU 103 (Mutato)"
puts "   - Viola: Ponte H inter-repeat protetto HIS 40 -> SER 74"
puts " Per salvare l'immagine per la Slide 6: File -> Render -> Snapshot / Tachyon"
puts "==========================================================================\n"
