"""Figure 1L - genome-track example PF3D7_1334800 (missing 5′ UTR defined by long reads)
Standalone: python panel_L.py   (reads ../data/Figure1_L_*.tsv; writes ../panels/Figure1_L.pdf/.svg/.png)
Panel letter and track row labels are NOT drawn (letter added in Illustrator; row labels = Figure1_JKL_track_labels)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from track_panel import draw_track_panel
draw_track_panel('L')
