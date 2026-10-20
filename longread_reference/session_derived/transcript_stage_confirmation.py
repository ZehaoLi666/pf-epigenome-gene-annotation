"""Session-derived computation of transcript_stage_confirmation.tsv.gz.

Recovered from the record of the interactive analysis session in which transcript_stage_confirmation.tsv.gz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: hc_reference.gtf.gz, stage_TPM_9stages.tsv, loci_transcript_classes.tsv, tss_clusters.tsv.gz, hc_model_evidence.tsv.gz, hc_class_min.tsv, junction_uniq_counts_by_stage.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import gzip
import re
import collections

ev = pd.read_csv(result('hc_model_evidence.tsv.gz'), sep="\t", low_memory=False)
tc = pd.read_csv(result('tss_clusters.tsv.gz'), sep="\t", low_memory=False)
sq = pd.read_csv(result('hc_class_min.tsv'), sep="\t")
jc = pd.read_csv(result('junction_uniq_counts_by_stage.tsv'), sep="\t")
jc = jc.set_index(["chrom", "start", "end"])
U = pd.read_csv(result('loci_transcript_classes.tsv'), sep="\t")
T = pd.read_csv(result('stage_TPM_9stages.tsv'), sep="\t", index_col=0)

ex = collections.defaultdict(list)
with gzip.open(result('hc_reference.gtf.gz'), "rt") as fh:
    for ln in fh:
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        ex[re.search('transcript_id "([^"]+)"', f[8]).group(1)].append((int(f[3]), int(f[4])))

SRS = {"ring": "Ring", "troph": "Trophozoite", "schizont": "Schizont"}
jd = {s: jc[SRS[s]].to_dict() for s in SRS}
jpool = jc.sum(axis=1).to_dict()

def jsup(tid, ch, s, minr=2):
    e = sorted(ex[tid])
    ints = [(ch, e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1)]
    if not ints:
        return True
    d = jd[s] if s else jpool
    return all(d.get(k, 0) >= minr for k in ints)

ev = ev.merge(sq.rename(columns={"isoform": "transcript_id"}), on="transcript_id", how="left")
ev["multi_exon"] = ev.n_exons > 1
for s in SRS:
    ev[f"sj_{s}"] = [jsup(t, ch, s) for t, ch in zip(ev.transcript_id, ev.chrom)]
ev["sj_any"] = [jsup(t, ch, None) for t, ch in zip(ev.transcript_id, ev.chrom)]

tcx = tc.assign(models=tc.models.astype(str).str.split(",")).explode("models")
for s in SRS:
    tcx[f"cage2_{s}"] = (tcx[f"cage_{s}_A"] >= 1) & (tcx[f"cage_{s}_B"] >= 1)
m2c = tcx.drop_duplicates("models").set_index("models")
for s in SRS:
    ev[f"cage_{s}_2rep"] = ev.transcript_id.map(m2c[f"cage2_{s}"]).fillna(False).astype(bool)

ev["struct_ok"] = ~ev.structural_category.isin(["fusion", "genic_intron"]) & (ev.RTS_stage.astype(str).str.upper() != "TRUE")
ev["tes_ok"] = (ev.tes_internal_priming_flag == 0) & (ev.tes_reads >= 2)
cls = U.set_index("locus").transcript_class
ev["unit_class"] = ev.locus.map(cls).fillna("other")
ev["cgroup"] = ev.unit_class.map(lambda c: "mRNA" if c == "mRNA" else ("lncRNA" if c.startswith(("lncRNA", "novel_ORF")) else "ncRNA/other"))

STG = ["Ring", "ET", "LT", "Trophozoite", "Schizont", "Gametocyte", "Gii", "GV", "Ookinete"]

for s in SRS:
    m1 = ev[f"detect_{s}"] == "detected"
    m2 = m1 & ev.struct_ok & ev[f"sj_{s}"]
    m3 = m2 & ev[f"cage_{s}_2rep"]
    m4 = m3 & ev.tes_ok
    ev[f"confirmed_{s}"] = m4
    gl = ev[m4].groupby("locus").size()
    tp = T[SRS[s]].reindex(gl.index).fillna(0)
    gc = gl[tp >= 1].index
    ev[f"gene_confirmed_{s}"] = ev.locus.isin(gc)

cols = ["transcript_id", "locus", "cgroup", "unit_class", "chrom", "start", "end", "strand", "n_exons", "tier", "structural_category", "struct_ok", "tes_ok",
        *[f"{p}_{s}" for s in SRS for p in ["FL", "detect", "sj", "cage", "confirmed", "gene_confirmed"] if f"{p}_{s}" in ev.columns or p == "cage"]]
cols = [c.replace("cage_ring", "cage_ring_2rep").replace("cage_troph", "cage_troph_2rep").replace("cage_schizont", "cage_schizont_2rep") if c.startswith("cage_") else c for c in cols]
cols = [c for c in dict.fromkeys(cols) if c in ev.columns]
ev[cols].to_csv("transcript_stage_confirmation.tsv.gz", sep="\t", index=False)
