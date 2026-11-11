# ==============================================================================
# viz_02_solvated_box.tcl
# Script grafico VMD per generare la Figura della Scatola d'Acqua (Fase 1 - Slide 3):
#   - Carica il sistema solvatato e ionizzato (2myo_wb_ionized.psf / .pdb)
#   - Sfondo bianco, proiezione ortografica, senza assi
#   - Proteina (Myotrophin) al centro colorata per struttura secondaria (NewCartoon)
#   - Ioni Na+ (Blu) e Cl- (Giallo/Verde) come sfere (VDW)
#   - Scatola d'acqua semitrasparente + Linee del dominio periodico (PBC Box)
#
# Esecuzione dalla cartella principale BML/:
#   vmd -e scripts/viz_02_solvated_box.tcl
# ==============================================================================

mol new 01_structures/2myo_wb_ionized.psf
mol addfile 01_structures/2myo_wb_ionized.pdb

color Display Background white
axes location Off
display projection Orthographic
display depthcue off
display antialias on
display ambientocclusion on

mol delrep 0 top

# 1. Proteina al centro in NewCartoon (colorata per Repeat / Index dal Blu al Rosso)
mol representation NewCartoon 0.32 20.0 4.1 0
mol color Index
mol selection "protein"
mol material AOChalky
mol addrep top

# 2. Ioni Sodio Na+ (SOD) -> Sfere Blu (ColorID 0)
mol representation VDW 0.8 20.0
mol color ColorID 0
mol selection "resname SOD"
mol material AOEdgy
mol addrep top

# 3. Ioni Cloro Cl- (CLA) -> Sfere Arancioni/Gialle (ColorID 3)
mol representation VDW 0.8 20.0
mol color ColorID 3
mol selection "resname CLA"
mol material AOEdgy
mol addrep top

# 4. Acqua (TIP3P) -> Reticolo sottile Ciano (Lines)
mol representation Lines 1.0
mol color ColorID 10
mol selection "water"
mol material Transparent
mol addrep top

# 5. Disegna il bordo geometrico della cella periodica (PBC Box)
package require pbctools
pbc set {71.840 58.972 54.840} -all
pbc box -center bb -color black -width 2.0

display resetview
puts "\n=================================================================="
puts "  SCENA 3D DELLA SCATOLA D'ACQUA PRONTA PER LA FOTO (Slide 3)!"
puts "  - Proteina al centro: NewCartoon (gradiente N-term -> C-term)"
puts "  - Sfere Blu (23):     Ioni Sodio (Na+)"
puts "  - Sfere Arancio (18): Ioni Cloro (Cl-)"
puts "  - Reticolo Ciano:     6535 molecole d'acqua TIP3P + Bordo PBC nero"
puts ""
puts "  Salva l'immagine con 'File -> Render -> Snapshot' in 06_presentation/"
puts "==================================================================\n"
