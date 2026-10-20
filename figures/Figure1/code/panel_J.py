"""Figure 1J - genome-track example PF3D7_0424600 (5′ UTR extends over an upstream lncRNA)
Standalone: python panel_J.py   (reads ../data/Figure1_J_*.tsv; writes ../panels/Figure1_J.pdf/.svg/.png)
Panel letter and track row labels are NOT drawn (letter added in Illustrator; row labels = Figure1_JKL_track_labels)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from track_panel import draw_track_panel
draw_track_panel('J')
