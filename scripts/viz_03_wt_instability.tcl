# ==============================================================================
# viz_03_wt_instability.tcl
# Fase 2 (ROAD2.md - Foto Slide 4):
# Visualizza in 3D la Myotrophin WT dopo stress termico a 400 K:
#   - Struttura iniziale nativa (Frame 0) in Grigio trasparente (riferimento)
#   - Struttura finale a 400 K colorata in base a <Q_i> (campo Beta):
#       Rosso = Residui che hanno perso i contatti nativi (<Q_i> basso, instabili)
#       Blu   = Residui stabili (<Q_i> vicino a 1.0)
#
# Esecuzione dalla cartella principale BML/:
#   vmd -e scripts/viz_03_wt_instability.tcl
# ==============================================================================

mol new 02_wt_baseline/wt_initial_ref.pdb
mol new 02_wt_baseline/wt_400K_q_beta.pdb

color Display Background white
axes location Off
display projection Orthographic
display depthcue off
display antialias on
display ambientocclusion on

# 1. Mol 0: Struttura iniziale (riferimento nativo) in Grigio semitrasparente
mol delrep 0 0
mol representation NewCartoon 0.26 20.0 4.1 0
mol color ColorID 2
mol selection "protein"
mol material Transparent
mol addrep 0

# 2. Mol 1: Struttura a 400 K colorata per <Q_i> (salvato nella colonna Beta, scala RWB: Rosso=0.4 .. Blu=1.0)
mol delrep 0 1
color scale method RWB
mol representation NewCartoon 0.32 20.0 4.1 0
mol color Beta
mol selection "protein"
mol material AOChalky
mol addrep 1
mol scaleminmax 1 0 0.45 0.95

# 3. Evidenzia in Licorice i residui non-mascherati piu' instabili (Beta < 0.65)
mol representation Licorice 0.22 20.0 20.0
mol color Beta
mol selection "protein and noh and beta < 0.65 and not resid 1 2 3 116 117 118 4 5 7 8 9 11 12 14 16 17 34 36 41 44 45 47 66 67 69 77 78 80 99 100"
mol material AOEdgy
mol addrep 1
mol scaleminmax 1 1 0.45 0.95

display resetview
puts "\n=================================================================="
puts "  SCENA 3D INSTABILITA' WT A 400 K PRONTA PER LA FOTO (Slide 4)!"
puts "  - Grigio Trasparente: Posizione nativa iniziale (Frame 0)"
puts "  - Colorato (Rosso -> Blu): Struttura a 400 K colorata per <Q_i>"
puts "      * ROSSO = Punti deboli che perdono contatti (<Q_i> basso)"
puts "      * BLU   = Core stabile (<Q_i> alto)"
puts "==================================================================\n"
