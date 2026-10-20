"""Normalisation, per-stage labels/masks and annotation helpers for the Pf boundary run."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'boundary_model')
from bm.tracks import RES
from bm import labels as L

import os
STAGES_REAL = ["ring", "troph", "schizont"]
STAGES = ["pool"] if os.environ.get("POOL") else STAGES_REAL
META = json.load(open("data/meta.json"))
CS = META["chrom_sizes"]; CHROMS = sorted(CS)
TEST = [c for c in CHROMS if 1 <= int(c[6:8]) <= 5]
VAL = [c for c in CHROMS if int(c[6:8]) in (6, 7)]
TRAIN = [c for c in CHROMS if int(c[6:8]) >= 8]
LAB = LAB_DIR
A = {  # result tables (long-read reference v1, confirmation tables)
 'tsc': 'transcript_stage_confirmation.tsv.gz', 'hce': 'hc_model_evidence.tsv.gz',
 'etc': 'end_transcription_contexts.tsv.gz', 'gsc': 'gene_stage_confirmation.tsv',
 'ltc': 'loci_transcript_classes.tsv', 'tpm': 'stage_TPM_9stages.tsv',
 'alt': 'confirmed_alternative_TSSs.tsv', 'sw': 'confirmed_genes_stage_TSS_switches.tsv'}


def load_raw(kind, stage):
    z = np.load(f"data/raw_{kind}_{stage}.npz"); return {c: z[c] for c in CHROMS}


def normalise_stage(Xs: dict, tracks: list, quantile_tracks=("H3K9me3",)):
    """log1p + z-score per (stage, track), stats from TRAIN chromosomes of that stage.
    Tracks in quantile_tracks are first rank-transformed genome-wide within stage
    (cross-study batch: H3K9me3 comes from three different papers)."""
    out = {}
    for s, X in Xs.items():
        Y = {c: np.log1p(np.clip(v, 0, None)) for c, v in X.items()}
        for t in quantile_tracks:
            if t not in tracks: continue
            i = tracks.index(t)
            allv = np.concatenate([Y[c][i] for c in CHROMS])
            srt = np.sort(allv)
            for c in CHROMS:
                Y[c][i] = np.searchsorted(srt, Y[c][i], side="right") / len(srt)
        tr = np.concatenate([Y[c] for c in TRAIN], 1)
        mu, sd = tr.mean(1, keepdims=True), tr.std(1, keepdims=True) + 1e-6
        out[s] = {c: ((Y[c] - mu) / sd).astype(np.float32) for c in CHROMS}
    return out


def cluster_ends(df, gap=50):
    """Merge ends of the same (chrom,strand,end_type) within `gap` bp -> one boundary
    (weighted median position by support w)."""
    rows = []
    for (c, st, et), g in df.sort_values("pos").groupby(["chrom", "strand", "end_type"]):
        p = g.pos.values; w = g.w.values; brk = np.r_[0, np.where(np.diff(p) > gap)[0] + 1, len(p)]
        for a, b in zip(brk[:-1], brk[1:]):
            pp, ww = p[a:b], w[a:b] + 1e-9
            o = np.argsort(pp); cw = np.cumsum(ww[o]); med = pp[o][np.searchsorted(cw, cw[-1] / 2)]
            rows.append(dict(chrom=c, strand=st, end_type=et, pos=int(med), n_models=b - a,
                             support=float(w[a:b].sum()), iqr_bp=float(max(np.subtract(*np.percentile(pp, [75, 25])), 20)),
                             genes=",".join(sorted(set(g.gene.values[a:b])))[:200],
                             klass=",".join(sorted(set(g.klass.values[a:b])))))
    return pd.DataFrame(rows)


def transcript_ends():
    hce = pd.read_csv(host_path('hce'), sep='\t', low_memory=False)
    tsc = pd.read_csv(host_path('tsc'), sep='\t', low_memory=False)
    d = hce.merge(tsc[["transcript_id", "unit_class"] + [f"confirmed_{s}" for s in STAGES]], on="transcript_id")
    d = d[d.chrom.isin(CHROMS)]
    rows = []
    for et, col in (("TSS", "tss"), ("TTS", "tes")):
        e = d[["transcript_id", "chrom", "strand", col, "locus", "unit_class", "n_reads"] + [f"confirmed_{s}" for s in STAGES]].rename(columns={col: "pos"})
        e["end_type"] = et; rows.append(e)
    return pd.concat(rows, ignore_index=True)


_HP = {}
def host_path(k):
    if k not in _HP:
        import host_paths; _HP.update(host_paths.P)
    return _HP[k]


def stage_reference(ends, stage):
    """gold = clustered ends of transcripts confirmed in `stage`; bronze = clustered ends of
    other high-confidence models (real-but-unconfirmed -> masked, never negative)."""
    e = ends.assign(w=ends.n_reads.clip(lower=1), gene=ends.locus.astype(str), klass=ends.unit_class.astype(str))
    conf = e[e[f"confirmed_{stage}"] == True]
    g = cluster_ends(conf); g["tier"] = "gold"
    b = cluster_ends(e[e[f"confirmed_{stage}"] != True]); b["tier"] = "bronze"
    r = pd.concat([g, b], ignore_index=True)
    r["boundary_id"] = [f"{stage}_{i}" for i in range(len(r))]; r["gene_id"] = r.genes
    return r


def exclusion_regions(pad_v48=300, pad_lnc=300):
    """Channel-specific unlabelled regions: v48 gene starts -> TSS channel, gene ends -> TTS
    channel (+/-300 bp, strand-aware); same for catalogued lncRNA starts/ends."""
    rows = []
    def add(chrom, lo, hi, strand, pad):
        if strand not in "+-": 
            for st in "+-": add(chrom, lo, hi, st, pad)
            return
        tss, tts = (lo, hi) if strand == "+" else (hi, lo)
        rows.append((chrom, max(0, tss - pad), tss + pad, strand, "TSS"))
        if hi is not None and hi != lo: rows.append((chrom, max(0, tts - pad), tts + pad, strand, "TTS"))
    for r in gff_genes().itertuples(): add(r.chrom, r.start, r.end, r.strand, pad_v48)
    for f in ["3 lncRNA datasets/Karine_2023.bed", "3 lncRNA datasets/Marcus_2022.bed", "3 lncRNA datasets/Pardis_2015.bed"]:
        b = pd.read_csv(f"{LAB}/{f}", sep="\t", header=None, usecols=[0, 1, 2, 5])
        for r in b.itertuples(index=False): add(r[0], r[1], r[2], r[3], pad_lnc)
    b = pd.read_csv(f"{LAB}/1774 identified lncRNA/unique_lncRNA_TSS_1774.bed", sep="\t", header=None, usecols=[0, 1, 5])
    for r in b.itertuples(index=False):
        for st in ([r[2]] if r[2] in "+-" else ["+", "-"]):
            rows.append((r[0], max(0, r[1] - pad_lnc), r[1] + pad_lnc, st, "TSS"))
    x = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "end_type"])
    return x[x.chrom.isin(CHROMS)].reset_index(drop=True)


def gff_genes(version=48):
    f = f"{LAB}/PlasmoDB/PlasmoDB-{version}/PlasmoDB-{version}_Pfalciparum3D7.gff"
    g = pd.read_csv(f, sep="\t", comment="#", header=None, low_memory=False,
                    names=["chrom", "src", "type", "start", "end", "score", "strand", "phase", "attr"])
    g = g[g.type.str.endswith("gene") & g.chrom.isin(CHROMS)].copy()
    g["gene_id"] = g.attr.str.extract(r"ID=([^;]+)")[0]
    g["biotype"] = g.type; g["desc"] = g.attr.str.extract(r"description=([^;]+)")[0]
    g["start"] -= 1  # 0-based
    return g[["chrom", "start", "end", "strand", "gene_id", "biotype", "desc"]].reset_index(drop=True)


def build_stage_labels(ref, excl):
    Y, M = L.build_labels(ref, CS, res=RES, tiers=("gold",), exclusion=None, ignore_radius_bp=100)
    for r in excl.itertuples(index=False):
        ch = L.channel_of(r.end_type, r.strand); lo, hi = r.start // RES, r.end // RES + 1
        seg = Y[r.chrom][ch, lo:hi]
        M[r.chrom][ch, lo:hi] = np.where(seg > 0.05, M[r.chrom][ch, lo:hi], 0.0)
    return Y, M
