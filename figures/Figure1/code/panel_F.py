"""Figure 1F - Comparison rules (schematic; no source data)
Standalone: python panel_F.py   (reads ../data/(none - schematic); writes ../panels/Figure1_F.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

fig, (axR,) = panel_figure('F')
axR.set_xlim(0, 100); axR.set_ylim(0, 23); axR.axis('off')
axR.text(0, 24.5, "Each v71 UTR end of the 2,389 reference genes was compared with every reference end cluster (main and alternative, ±100 bp); 5′ end shown", fontsize=6.3, va='bottom')
cases = [("Agrees", 'agree', 3, [(3.2, 14), (5.0, 14), (3.1, 13.8)], None), ("Over-long", 'ol_none', 1, [(4.0, 14), (5.2, 14), (4.1, 13.8)], (1, 4.0)),
         ("Too short", 'short', 4.5, [(2.5, 14), (2.8, 14), (2.6, 13.8)], (2.5, 4.5)), ("No annotated UTR", 'noutr', 6, [(3.5, 14), (3.8, 14), (3.6, 13.8)], (3.5, 6)),
         ("Between ends", 'between', 3.6, [(2.2, 14), (5.0, 14), (2.4, 13.8)], None)]
for i, (nm, ck, u5, reads, shade) in enumerate(cases):
    x0 = 9 + i * 18.4
    axR.text(x0 + 8.5, 19.5, nm, fontsize=6.3, fontweight='bold', ha='center', color=COL[ck] if ck != 'ol_none' else '#c0640f')
    yy = 13.0; axR.plot([x0, x0 + 17], [yy, yy], color='#bbbbbb', lw=0.5)
    if u5 < 6: axR.add_patch(Rectangle((x0 + u5, yy - 1.1), 6 - u5, 2.2, fc=COL['utr'], ec='none'))
    axR.add_patch(Rectangle((x0 + 6, yy - 2.0), 6, 4.0, fc=COL['cds'], ec='none')); axR.add_patch(Rectangle((x0 + 12, yy - 1.1), 2, 2.2, fc=COL['utr'], ec='none'))
    for kk, (p, q) in enumerate(reads): axR.plot([x0 + p, x0 + q], [7.5 - kk * 1.9] * 2, color='#8c8c8c', lw=1.0, solid_capstyle='butt')
    if shade: axR.add_patch(Rectangle((x0 + shade[0], 2.0), shade[1] - shade[0], 13.5, fc=COL[ck], alpha=0.4, ec='none', zorder=0))
axR.text(0.5, 13.0, "PlasmoDB v71", fontsize=5.5, va='center'); axR.text(0.5, 5.6, "reference\nend clusters", fontsize=5.5, va='center', linespacing=0.9)
save_panel(fig, 'F')
