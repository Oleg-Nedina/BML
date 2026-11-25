# ==============================================================================
# viz_01_structure_and_mask.tcl
# Script grafico VMD per generare le Figure della Slide 2 (Fase 1 - ROAD1.md):
#   - Sfondo bianco pulito da pubblicazione, proiezione ortografica, senza assi
#   - I 4 Ankyrin Repeats (ANK1..ANK4) colorati con 4 colori distinti in NewCartoon
#   - La triade funzionale chiave (TRP 8, ASP 44, GLU 78) in VDW/Licorice Rosso
#   - Le 2 Istidine chiave inter-repeat (HIS 40 e HIS 67) e i loro partner H-bond
#     (SER 74, LYS 37, GLY 35) in Licorice Ciano/Giallo
#
# Esecuzione dalla cartella principale BML/:
#   vmd -e scripts/viz_01_structure_and_mask.tcl
# ==============================================================================

mol new 01_structures/raw/2MYO.pdb

# 1. Impostazioni grafiche da pubblicazione (sfondo bianco, ortografico, luce pulita)
color Display Background white
axes location Off
display projection Orthographic
display depthcue off
display antialias on
display ambientocclusion on

# Cancella la rappresentazione di default (Lines)
mol delrep 0 top

# 2. Rappresentazione 0: N-terminal (resid 1 to 27: H1, H2) - Blu (ColorID 0)
mol representation NewCartoon 0.30 20.0 4.1 0
mol color ColorID 0
mol selection "protein and resid 1 to 27"
mol material AOChalky
mol addrep top

# 3. Rappresentazione 1: Ankyrin Repeat 1 (resid 28 to 60: loop 32-36, H3, H4) - Ciano/Teal (ColorID 10)
mol representation NewCartoon 0.30 20.0 4.1 0
mol color ColorID 10
mol selection "protein and resid 28 to 60"
mol material AOChalky
mol addrep top

# 4. Rappresentazione 2: Ankyrin Repeat 2 (resid 61 to 93: loop 65-69, H5, H6) - Arancione (ColorID 3)
mol representation NewCartoon 0.30 20.0 4.1 0
mol color ColorID 3
mol selection "protein and resid 61 to 93"
mol material AOChalky
mol addrep top

# 5. Rappresentazione 3: C-terminal (resid 94 to 118: H7) - Verde (ColorID 7)
mol representation NewCartoon 0.30 20.0 4.1 0
mol color ColorID 7
mol selection "protein and resid 94 to 118"
mol material AOChalky
mol addrep top

# 6. Rappresentazione 4: Triade Funzionale Capping-Protein (TRP 8, ASP 44, GLU 78) - Rosso (ColorID 1)
mol representation Licorice 0.25 20.0 20.0
mol color ColorID 1
mol selection "protein and noh and resid 8 44 78"
mol material AOEdgy
mol addrep top

# 7. Rappresentazione 5: Tutti i 24 residui della Functional Mask (superficie trasparente rossa)
mol representation Surf 1.4 0.0
mol color ColorID 1
mol selection "protein and noh and resid 4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100"
mol material Transparent
mol addrep top
# Nascondiamo di default la superficie trasparente (puoi attivarla con doppio click in Graphics -> Representations!)
mol showrep top 5 off

# 8. Rappresentazione 6: Rete H-bond di HIS 40 (con SER 74) e HIS 67 (con LYS 37 e GLY 35) per amino-acid coloring
mol representation Licorice 0.20 20.0 20.0
mol color Name
mol selection "protein and noh and resid 35 37 40 67 74"
mol material AOShiny
mol addrep top

display resetview
puts "\n=================================================================="
puts "  SCENA 3D PRONTA PER LE FOTO DELLE SLIDE (Fase 1 - Yang et al. 1998)!"
puts "  - N-terminal       (1-27,  H1-H2): Blu"
puts "  - Ankyrin Repeat 1 (28-60, H3-H4): Ciano"
puts "  - Ankyrin Repeat 2 (61-93, H5-H6): Arancione"
puts "  - C-terminal       (94-118, H7):   Verde"
puts "  - Triade Funzionale (TRP 8, ASP 44, GLU 78): Bastoncini Rossi"
puts "  - Rete H-bond Istidine (HIS 40-SER 74 e LYS 37-HIS 67-GLY 35): Colori Atomici"
puts ""
puts "  Suggerimento: Apri 'Graphics -> Representations' per accendere/spegnere"
puts "  con un doppio click i singoli livelli e fare screenshot diversi!"
puts "  Per salvare: 'File -> Render -> Snapshot' (o Tachyon) in 06_presentation/"
puts "==================================================================\n"
