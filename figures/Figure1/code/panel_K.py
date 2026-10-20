"""Figure 1K - genome-track example PF3D7_0627500 (3′ UTR extends past the long-read 3′ ends)
Standalone: python panel_K.py   (reads ../data/Figure1_K_*.tsv; writes ../panels/Figure1_K.pdf/.svg/.png)
Panel letter and track row labels are NOT drawn (letter added in Illustrator; row labels = Figure1_JKL_track_labels)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from track_panel import draw_track_panel
draw_track_panel('K')
