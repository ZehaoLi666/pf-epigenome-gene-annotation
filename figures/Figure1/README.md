# Figure 1

Long-read reference and comparison with PlasmoDB v71 ends. py). `code/original/` keeps the scripts that drew the submitted composite (make_fig1.py, tracks.py, common.py); prep_data.py copies their derivations.

## Run

```bash
cd figures/Figure1
python code/prep_data.py   # inputs -> data/Figure1_*.tsv (needs pyBigWig and config/paths.yaml)
python code/run_all.py     # data/ -> panels/*.pdf|svg|png (numpy, pandas, matplotlib, scipy)
```

The per-panel source tables in `data/` are deposited with the paper; with them, `run_all.py` runs without any
other input. Panels are drawn at final size without panel letters; the figure is assembled in Adobe Illustrator.

## Code

| File | Role |
|---|---|
| `panel_A.py` | A: Long-read datasets and reference construction |
| `panel_B.py` | B: Read-level long-read coverage |
| `panel_C.py` | C: Completeness of stage-confirmed mRNA models |
| `panel_D.py` | D: Non-templated 5′ G at TSS clusters |
| `panel_E.py` | E: Internal-priming signature at 3′ ends |
| `panel_F.py` | F: Comparison rules (schematic) |
| `panel_G.py` | G: 5′ × 3′ class matrix |
| `panel_H.py` | H: 5′ UTR length vs farthest reference end |
| `panel_I.py` | I: 3′ UTR length vs farthest reference end |
| `panel_J.py` | J: Example PF3D7_0424600 |
| `panel_JKL_track_labels.py` | draws panel JKL track labels |
| `panel_K.py` | K: Example PF3D7_0627500 |
| `panel_L.py` | L: Example PF3D7_1334800 |
| `prep_data.py` | gathers / derives all inputs, writes data/ |
| `run_all.py` | runs every panel script |
| `style.py` | shared rcParams, colours, panel geometry |
| `track_panel.py` | shared genome-track drawing routine (J–L) |
| `original/` | scripts that drew the submitted composite (provenance) |

## Inputs

Intermediate result tables (`RESULTS_DIR`):

- `v71_vs_reference_ends.tsv`
- `funnel.tsv`
- `genome_coverage.tsv`
- `gene_longread_coverage.tsv.gz`
- `model_completeness.tsv`
- `tss_read_evidence_by_group.tsv`
- `internal_priming_v71.tsv`
- `gene_boundaries_v3_all_genes.tsv`

Track / annotation files (`LAB_DIR`):

- `PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff`
- `IGV_longread_transcript_models/confirmed_transcripts_by_category.bed`
- `IGV_longread_transcript_models/lncRNA_catalogue_3datasets.bed`
- `IGV_longread_transcript_models/DAFTseq_pool_{plus,minus}_CPM_10bp.bw`
- `lncRNA project tracks/cage_tss/CTSS_<stage>_{pos,neg}_CPM.bw`

Checked for this release: `prep_data.py` reproduced every file in `data/` byte-identically and `run_all.py` reproduced every PNG panel pixel-identically.
