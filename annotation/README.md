# Refined annotation

`build_resource.py` builds `refined_utr_annotation.gff3` / `.tsv` (Data S1) from the per-gene boundary table
(`gene_boundaries_v3_all_genes.tsv`; long-read ends where available, otherwise the model) and the PlasmoDB-68
GFF (CDS / exon structure; identical to v71). `validate_resource.py` checks the GFF3 (feature counts, coordinate
order, chromosome bounds). Both `exec` `tracks.py` (GFF / BED readers) and `common.py` (plot style) from this folder:
run them from a folder where `fig6/tracks.py` and `fig6/common.py` point to these files, or edit the `exec` lines.
