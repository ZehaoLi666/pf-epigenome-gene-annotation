# Figure 3

Boundary model and cross-validated performance against long-read ends.

## Run

```bash
cd figures/Figure3
python code/prep_data.py   # inputs -> data/Figure3_*.tsv (needs pyBigWig and config/paths.yaml)
python code/run_all.py     # data/ -> panels/*.pdf|svg|png (numpy, pandas, matplotlib, scipy)
```

The per-panel source tables in `data/` are deposited with the paper; with them, `run_all.py` runs without any
other input. Panels are drawn at final size without panel letters; the figure is assembled in Adobe Illustrator.

## Code

| File | Role |
|---|---|
| `legend_key.py` | legend key panel |
| `panel_A.py` | A: Boundary-model schematic |
| `panel_B.py` | B: Three-fold chromosome cross-validation |
| `panel_C.py` | C: Three example loci |
| `panel_D.py` | D: Cumulative distance to the nearest long-read end |
| `panel_E.py` | E: Input ablation vs long-read ends |
| `panel_F.py` | F: Full vs no-chromatin model by DAFT-seq expression quintile |
| `panel_G.py` | G: Where predictions fall (±200 bp) |
| `prep_data.py` | gathers / derives all inputs, writes data/ |
| `run_all.py` | runs every panel script |
| `style.py` | shared rcParams, colours, panel geometry |

## Inputs

Intermediate result tables (`RESULTS_DIR`):

- `gene_level_eval.tsv`
- `example_confirmed_models.tsv`
- `fig3E_ablation_longread.tsv`
- `fig3F_chromatin_by_expression_longread.tsv`

Track / annotation files (`LAB_DIR`):

- `PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff`
- `IGV_longread_transcript_models/DAFTseq_pool_{plus,minus}_CPM_10bp.bw`
- `lncRNA project tracks/cage_tss/CTSS_<stage>_{pos,neg}_CPM.bw`

Checked for this release: `prep_data.py` reproduced every file in `data/` byte-identically and `run_all.py` reproduced every PNG panel pixel-identically.
