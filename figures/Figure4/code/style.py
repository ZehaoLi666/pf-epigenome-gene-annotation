"""Shared style for the Figure 4 panel scripts.

Fonts, colours and line widths are copied from the code that produced the approved composite
(make_fig4.py + common.py). Every panel is drawn on its own
canvas whose size equals the panel's bounding box in the 7.2 x 8.9 in composite, and every axes keeps
its original position in inches, so each panel file drops into Illustrator at 100 % scale.
Panel letters (A-F) are NOT drawn; their composite positions are listed in LETTER_POS_MM / README.md.
"""
from pathlib import Path
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'data'
PANELS = HERE.parent / 'panels'

# ---- rcParams -------------------------------------------------------------------------------------
# 1) defaults of the lab's figure-style helper (apply_figure_style) ...
FIGURE_STYLE = {'axes.labelsize': 8.0, 'axes.linewidth': 0.6, 'axes.spines.right': False, 'axes.spines.top': False,
                'axes.titlelocation': 'left', 'axes.titlesize': 8.0, 'figure.dpi': 200.0, 'font.size': 8.0,
                'legend.fontsize': 7.0, 'legend.frameon': False, 'lines.linewidth': 1.2, 'patch.linewidth': 0.6,
                'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.dpi': 300.0,
                'xtick.labelsize': 6.0, 'xtick.major.size': 3.0, 'xtick.major.width': 0.6, 'ytick.labelsize': 6.0,
                'ytick.major.size': 3.0, 'ytick.major.width': 0.6}
# 2) ... overridden by the Figure 4 settings from common.py (these are what the composite uses) ...
COMMON = {'font.family': 'Arial', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5, 'xtick.labelsize': 5.5,
          'ytick.labelsize': 5.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
          'xtick.major.size': 2, 'ytick.major.size': 2, 'axes.spines.top': False, 'axes.spines.right': False,
          'pdf.fonttype': 42}
# 3) ... plus export settings: editable text in PDF/SVG, fixed canvas (no tight bbox, the canvas IS the panel),
#    patch.linewidth restored to the matplotlib default 1.0 that the approved composite was rendered with
#    (legend swatches in panel D), deterministic file output.
EXPORT = {'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.bbox': 'standard',
          'savefig.dpi': 300, 'patch.linewidth': 1.0, 'svg.hashsalt': 'figure4', 'figure.dpi': 100}
mpl.rcParams.update(FIGURE_STYLE); mpl.rcParams.update(COMMON); mpl.rcParams.update(EXPORT)

# ---- colours (as in the approved composite) ---------------------------------------------------------
COL = dict(agree='#4a78c2', ol_unit='#a6400a', ol_lnc='#e08214', ol_anti='#f5b041', ol_none='#fbd9a6', short='#7b4fa0',
           noutr='#2a9d8f', between='#c4a46b', na='#e3e3e3', ref='#1f3b73', cds='#1f3b73', utr='#8fb3e0', lnc='#e08214',
           cage='#c23b6e', daft='#4a78c2', tss='#7b4fa0', tts='#2a9d8f', body='#4a78c2', inter='#b5b5b5',
           trunc='#9a9a9a', art='#d9822b')
REF, V71 = COL['ref'], COL['art']          # refined annotation (navy) / PlasmoDB v71 (orange)
EXPT = '#4d4d4d'                             # long-read reference (experimental)
LIGHT = '#7fa7d9'                            # model without chromatin
COLM = {'Ours_chromatin': REF, 'Ours_nochromatin': LIGHT, 'FixedOffset': '#9a9a9a', 'v71all': V71,
        'ChromHMM': '#7fb685', 'StringTie_short': '#8c6d46', 'StringTie_long': '#b8955a'}
SM = 5.3                                     # small annotation font size (pt)
TITLE_FS, LETTER_FS = 7, 9

# ---- composite geometry ------------------------------------------------------------------------------
FIG_W, FIG_H, DPI = 7.2, 8.9, 300            # composite size (in) and raster resolution
# Panel boxes in composite pixels at 300 dpi (x0, y0, x1, y1; origin bottom-left) = union of the panel's axes,
# tick labels, legends, title and letter, padded by 6 px (0.5 mm).
REGIONS_PX = {'A': (149, 2015, 1016, 2621), 'B': (1261, 2015, 2128, 2621), 'C': (52, 1417, 2128, 1985),
              'D': (91, 776, 1063, 1330), 'E': (1289, 774, 2123, 1330), 'F': (80, 10, 2112, 613)}
# Panel titles and letters: (x, y) in composite figure fractions, ha='left', va='bottom'.
TITLES = {'A': ((0.09672376543209876, 0.9675000000000001), 'TSS accuracy (5′-end consensus)'),
          'B': ((0.6117237654320987, 0.9675000000000001), 'TTS accuracy (Siegel 2014)'),
          'C': ((0.05207445008077739, 0.7295), 'Independent 5′-end and poly(A) data centred on refined vs v71 ends (% within ±100 bp at right)'),
          'D': ((0.0693837962962963, 0.4838250343808521), 'Refined ends relative to v71'),
          'E': ((0.623783950617284, 0.4838250343808521), 'UTR lengths'),
          'F': ((0.0642854938271605, 0.21543033707865167), 'Accuracy by expression level: reference genes and genes whose ends come from the model')}
LETTERS = {'A': (0.07272376543209877, 0.9660000000000002), 'B': (0.5877237654320987, 0.9660000000000002),
           'C': (0.028074450080777396, 0.7280000000000001), 'D': (0.045383796296296304, 0.4823250343808521),
           'E': (0.599783950617284, 0.4823250343808521), 'F': (0.0402854938271605, 0.21393033707865167)}


def region_mm(letter):
    x0, y0, x1, y1 = REGIONS_PX[letter]
    return (x1 - x0) / DPI * 25.4, (y1 - y0) / DPI * 25.4


def letter_pos_mm(letter):
    """Panel-letter baseline-left position in mm from the panel file's TOP-LEFT corner (9 pt bold Arial)."""
    x0, y0, x1, y1 = REGIONS_PX[letter]
    lx, ly = LETTERS[letter]
    return (lx * FIG_W * DPI - x0) / DPI * 25.4, (y1 - ly * FIG_H * DPI) / DPI * 25.4


class Panel:
    """A canvas the size of one composite panel; composite figure-fraction coordinates are mapped onto it."""

    def __init__(self, letter, title=True):
        self.letter = letter
        x0, y0, x1, y1 = REGIONS_PX[letter]
        self.x0, self.y0 = x0 / DPI, y0 / DPI               # panel origin in composite inches
        self.w, self.h = (x1 - x0) / DPI, (y1 - y0) / DPI    # panel size in inches
        self.fig = plt.figure(figsize=(self.w, self.h))
        if title:
            (tx, ty), ttl = TITLES[letter]
            self.fig.text(*self.xy(tx, ty), ttl, fontsize=TITLE_FS, ha='left', va='bottom')

    def xy(self, x, y):
        return ((x * FIG_W - self.x0) / self.w, (y * FIG_H - self.y0) / self.h)

    def rect(self, r):
        l, b, w, h = r
        return [(l * FIG_W - self.x0) / self.w, (b * FIG_H - self.y0) / self.h, w * FIG_W / self.w, h * FIG_H / self.h]

    def add_axes(self, r):
        return self.fig.add_axes(self.rect(r))

    def save(self, stem=None):
        stem = stem or f'Figure4_{self.letter}'
        PANELS.mkdir(exist_ok=True)
        self.fig.savefig(PANELS / f'{stem}.pdf', metadata={'CreationDate': None, 'ModDate': None})
        self.fig.savefig(PANELS / f'{stem}.svg', metadata={'Date': None})
        self.fig.savefig(PANELS / f'{stem}.png', dpi=DPI)
        plt.close(self.fig)
        w, h = self.w * 25.4, self.h * 25.4
        print(f'{stem}: {w:.1f} x {h:.1f} mm -> panels/{stem}.pdf|.svg|.png')
