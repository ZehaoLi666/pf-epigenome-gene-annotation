# Figure 5

Promoter architecture, head-to-head gene pairs and stage-specific TSS usage. Statistics behind printed p-values are written to data/Figure5_statistics.tsv.

## Run

```bash
cd figures/Figure5
python code/prep_data.py   # inputs -> data/Figure5_*.tsv (needs pyBigWig and config/paths.yaml)
python code/run_all.py     # data/ -> panels/*.pdf|svg|png (numpy, pandas, matplotlib, scipy)
```

The per-panel source tables in `data/` are deposited with the paper; with them, `run_all.py` runs without any
other input. Panels are drawn at final size without panel letters; the figure is assembled in Adobe Illustrator.

## Code

| File | Role |
|---|---|
| `legend_AD.py` | legend for panels A–D |
| `panel_A.py` | A: ATAC metagene at TSS |
| `panel_B.py` | B: ApiAP2 metagene at TSS |
| `panel_C.py` | C: MNase metagene at TSS |
| `panel_D.py` | D: GTGCA motif density |
| `panel_E.py` | E: Neighbouring-gene co-expression by orientation |
| `panel_F.py` | F: MNase signal across head-to-head pairs |
| `panel_G.py` | G: Head-to-head pairs by distance and NDR |
| `panel_H.py` | H: Stage-specific TSS switching: Grade A/B pairs vs single-TSS controls |
| `panel_I.py` | I: Focal genes |
| `panel_J.py` | J: uAUG counts in Grade A/B TSS pairs |
| `prep_data.py` | gathers / derives all inputs, writes data/ |
| `run_all.py` | runs every panel script |
| `style.py` | shared rcParams, colours, panel geometry |

## Inputs

Intermediate result tables (`RESULTS_DIR`):

- `gene_boundaries_v3_all_genes.tsv`
- `gene_daft_expression_truth.tsv.gz`
- `apiap2_genebody_track_summary.tsv`
- `bidirectional_pairs.tsv`
- `divergent_pairs_mnase_ndr.tsv`
- `tss_decomposition_by_pair.tsv`
- `tss_decomposition_controls.tsv`
- `MNase-seq_T{5,10,...,40}_qNorm.bw (Kensche 2016 MNase-seq, nucDetective qNorm tracks, Zenodo 16779899)`

Track / annotation files (`LAB_DIR`):

- `IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta`
- `lncRNA project tracks/accessibility_nucleosome/ (ATAC)`
- `lncRNA project tracks/apiap2_tf/ (ApiAP2 ChIP-seq)`

Checked for this release: `prep_data.py` reproduced every file in `data/` byte-identically and `run_all.py` reproduced every PNG panel pixel-identically.
