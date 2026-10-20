IGV view of long-read confirmed genes and transcript models (PlasmoDB v68 coordinates; v48 = v68 nuclear sequence)
Open in IGV: File > Open Session... > igv_session.xml (keep all files in this folder together).

Tracks
- PlasmoDB v68 genes: reference annotation (GFF3, from your PlasmoDB-68 folder).
- LR confirmed genes: 2,776 v68 protein-coding genes confirmed by long reads (gene-level, >=1 stage); name = gene|n models|full-CDS|5'trunc|3'trunc|lncRNA sense (confirmed/catalogue)|antisense.
  colour: dark blue single model = v68; green full-CDS end variants only; light blue full-CDS + truncated; purple full-CDS + novel splice; grey only truncated models.
- LR confirmed transcripts: 13,348 transcript models confirmed in >=1 stage (12,934 mRNA). name = id|category|locus|stages (R/T/S)|tier. colour:
  dark blue = full CDS, ends match v68 (+-200 bp); light blue = full CDS, TSS inside annotated 5'UTR; green = full CDS, TSS upstream of annotation;
  teal = full CDS, alternative 3' end only; purple = full CDS, novel splice structure; grey = starts inside CDS (5'-truncated);
  light grey = ends inside CDS (3'-truncated); magenta = intergenic lncRNA; orange = antisense lncRNA; brown = other ncRNA; black = mRNA locus absent from v68.
- lncRNA catalogues: the three published lncRNA BED files (name = source|id).
- Over-long 5'UTR: part of the annotated 5'UTR beyond the farthest long-read-supported TSS (469 genes); red = held-out ML model agrees (TSS at long-read end, excess intergenic), pink = other.
- CAGE: stage CTSS bigWigs (+/- strand), from your lncRNA project tracks folder.
