# Figure 2

Training labels.

## Run

```bash
cd figures/Figure2
python code/prep_data.py   # inputs -> data/Figure2_*.tsv (needs pyBigWig and config/paths.yaml)
python code/run_all.py     # data/ -> panels/*.pdf|svg|png (numpy, pandas, matplotlib, scipy)
```

The per-panel source tables in `data/` are deposited with the paper; with them, `run_all.py` runs without any
other input. Panels are drawn at final size without panel letters; the figure is assembled in Adobe Illustrator.

## Code

| File | Role |
|---|---|
| `panel_A.py` | A: Label-definition schematic |
| `panel_B.py` | B: Number of labels per class |
| `panel_C.py` | C: Input-track metagene coverage |
| `prep_data.py` | gathers / derives all inputs, writes data/ |
| `run_all.py` | runs every panel script |
| `style.py` | shared rcParams, colours, panel geometry |

## Inputs

Intermediate result tables (`RESULTS_DIR`):

- `ml3_end_positives.tsv`
- `flank_gaps_v3_longread_tiers.tsv`
- `metagene_profiles.npz`

Checked for this release: `prep_data.py` reproduced every file in `data/` byte-identically and `run_all.py` reproduced every PNG panel pixel-identically.
