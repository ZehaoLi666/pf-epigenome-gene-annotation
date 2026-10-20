# Figure 4

Benchmark against independent 5′- and 3′-end data and comparison with PlasmoDB v71. If the Figure4_*.tsv tables of the original composite are present in RESULTS_DIR, prep_data.py asserts that it reproduces them.

## Run

```bash
cd figures/Figure4
python code/prep_data.py   # inputs -> data/Figure4_*.tsv (needs pyBigWig and config/paths.yaml)
python code/run_all.py     # data/ -> panels/*.pdf|svg|png (numpy, pandas, matplotlib, scipy)
```

The per-panel source tables in `data/` are deposited with the paper; with them, `run_all.py` runs without any
other input. Panels are drawn at final size without panel letters; the figure is assembled in Adobe Illustrator.

## Code

| File | Role |
|---|---|
| `panel_A.py` | A: TSS accuracy vs 5′-end consensus |
| `panel_B.py` | B: TTS accuracy vs Siegel 2014 major poly(A) site |
| `panel_C.py` | C: 5′-end / poly(A) signal centred on refined vs v71 ends |
| `panel_D.py` | D: Refined ends relative to v71, by source |
| `panel_E.py` | E: 5′ / 3′ UTR length distributions |
| `panel_F.py` | F: Accuracy by DAFT-seq expression quintile |
| `prep_data.py` | gathers / derives all inputs, writes data/ |
| `run_all.py` | runs every panel script |
| `style.py` | shared rcParams, colours, panel geometry |

## Inputs

Intermediate result tables (`RESULTS_DIR`):

- `benchmark_method_predictions.tsv`
- `benchmark_truth_sets.tsv`
- `gene_boundaries_v3_all_genes.tsv`
- `siegel2014_pas_per_gene.tsv`
- `chappell2020_5UTRseq_tagclusters.tsv`
- `shaw2021_5CAPture_5CN_clusters.tsv`
- `daft_stratified_reference_vs_model.tsv`

Track / annotation files (`LAB_DIR`):

- `lncRNA project tracks/cage_tss/CTSS_<stage>_{pos,neg}_CPM.bw`

Checked for this release: `prep_data.py` reproduced every file in `data/` byte-identically and `run_all.py` reproduced every PNG panel pixel-identically.
