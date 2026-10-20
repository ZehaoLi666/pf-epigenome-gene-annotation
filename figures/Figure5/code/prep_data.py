"""Figure 5 - derive every per-panel source table from the primary inputs.

Writes ../data/Figure5_*.tsv (+ one .npz matrix for panel F) and the panel
geometry table. Panel scripts read only these files.

INPUTS below are intermediate result tables read from RESULTS_DIR (config/paths.yaml) and
track files under LAB_DIR (see README, section "Inputs").
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import os, re
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import wilcoxon, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data")
os.makedirs(OUT, exist_ok=True)

TRACKS = f"{LAB_DIR}/lncRNA project tracks"
IGV = f"{LAB_DIR}/IGV_longread_transcript_models"
A = result

def art(vid, local):
    return A(vid) if A else local

INPUTS = {
    # gene boundaries: refined (long-read / chromatin-free model) and v71 TSS, ATG
    "genes": art("gene_boundaries_v3_all_genes.tsv", "inputs/gene_boundaries_v3_all_genes.tsv"),
    # truth table incl. Consensus_5end__TSS and DAFT-seq time courses
    "truth": art("gene_daft_expression_truth.tsv.gz", "inputs/gene_daft_expression_truth.tsv.gz"),
    # ApiAP2 track classification (promoter-type vs gene-body)
    "ap2_class": art("apiap2_genebody_track_summary.tsv", "inputs/apiap2_genebody_track_summary.tsv"),
    # all adjacent gene pairs with orientation and DAFT-seq correlation
    "pairs": art("bidirectional_pairs.tsv", "inputs/bidirectional_pairs.tsv"),
    # head-to-head pairs with exact refined anchors and MNase NDR calls
    "pairs_mnase": art("divergent_pairs_mnase_ndr.tsv", "inputs/divergent_pairs_mnase_ndr.tsv"),
    # stage-specific TSS decomposition (65 Grade A/B pairs) and single-TSS controls
    "decomp": art("tss_decomposition_by_pair.tsv", "inputs/tss_decomposition_by_pair.tsv"),
    "decomp_ctrl": art("tss_decomposition_controls.tsv", "inputs/tss_decomposition_controls.tsv"),
    "fasta": f"{IGV}/PlasmoDB-68_Pfalciparum3D7_Genome.fasta",
    "atac_dir": f"{TRACKS}/accessibility_nucleosome",
    "ap2_dir": f"{TRACKS}/apiap2_tf",
}
# Kensche et al. 2016 paired-end MNase-seq, nucDetective qNorm tracks (Zenodo 16779899)
MNASE = {tp: art(aid, f"inputs/MNase-seq_{tp}_qNorm.bw") for tp, aid in {
    "T5": "MNase-seq_T5_qNorm.bw", "T10": "MNase-seq_T10_qNorm.bw",
    "T15": "MNase-seq_T15_qNorm.bw", "T20": "MNase-seq_T20_qNorm.bw",
    "T25": "MNase-seq_T25_qNorm.bw", "T30": "MNase-seq_T30_qNorm.bw",
    "T35": "MNase-seq_T35_qNorm.bw", "T40": "MNase-seq_T40_qNorm.bw"}.items()}

W1, W2 = 1000, 500          # metagene window: -1000..+499 bp around the TSS
WH = 1500                   # panel F half-window around the pair midpoint


def smooth_track(v, w=21):
    cs = np.concatenate([[0], np.cumsum(v, dtype=np.float64)]); h = w // 2; n = len(v)
    lo = np.clip(np.arange(n) - h, 0, n); hi = np.clip(np.arange(n) + h + 1, 0, n)
    return ((cs[hi] - cs[lo]) / (hi - lo)).astype(np.float32)


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


# ---------------------------------------------------------------- signals
b = pyBigWig.open(MNASE["T5"]); CH = dict(b.chroms()); b.close()
raw = {c: np.zeros(L, np.float32) for c, L in CH.items()}
for tp, p in MNASE.items():
    b = pyBigWig.open(p)
    for c in CH:
        raw[c] += np.nan_to_num(np.array(b.values(c, 0, CH[c]), dtype=np.float32)) / len(MNASE)
    b.close()
MN = {c: smooth_track(raw[c]) for c in CH}; del raw
MED = float(np.median(np.concatenate(list(MN.values()))))

ATAC = {c: np.zeros(CH[c], np.float32) for c in CH}
for st in ["Ring", "Trophozoite", "Schizont"]:
    b = pyBigWig.open(f"{INPUTS['atac_dir']}/ATAC_{st}_minusInput.bw"); ch = b.chroms()
    for c in CH:
        if c in ch:
            ATAC[c] += np.nan_to_num(np.array(b.values(c, 0, CH[c]), dtype=np.float32)) / 3
    b.close()

AS = pd.read_csv(INPUTS["ap2_class"], sep="\t")
PROM = AS[AS.track_class.str.contains("promoter", case=False)].track.tolist()
assert len(PROM) == 10, PROM
AP2 = {c: np.zeros(CH[c], np.float32) for c in CH}
for t in PROM:
    b = pyBigWig.open(f"{INPUTS['ap2_dir']}/{t}_minusInput.bw"); ch = b.chroms()
    for c in CH:
        if c in ch:
            AP2[c] += np.nan_to_num(np.array(b.values(c, 0, CH[c]), dtype=np.float32)) / len(PROM)
    b.close()
GM_ATAC = float(np.mean(np.concatenate([ATAC[c] for c in CH])))
GM_AP2 = float(np.mean(np.concatenate([AP2[c] for c in CH])))

SEQ = {}; nm = None; buf = []
for l in open(INPUTS["fasta"]):
    if l[0] == ">":
        if nm: SEQ[nm] = "".join(buf)
        nm = l[1:].split()[0].split("|")[0].strip(); buf = []
    else:
        buf.append(l.strip().upper())
SEQ[nm] = "".join(buf)

# ---------------------------------------------------------------- panels A-D
G = pd.read_csv(INPUTS["genes"], sep="\t")
TR = pd.read_csv(INPUTS["truth"], sep="\t")
G = G.merge(TR[["gene_id", "Consensus_5end__TSS"]], on="gene_id", how="left")
G["ref"] = G.lr_TSS.where(G.lr_TSS.notna(), G.chromfree_model_TSS)
G["sg"] = np.where(G.strand == "+", 1, -1)
G["dpos"] = (G.ref - G.annot_tss) * G.sg
trim = G[(G.dpos > 300) & G.ref.notna() & G.annot_tss.notna()].reset_index(drop=True)
assert len(trim) == 1236
trim[["gene_id", "chrom", "strand", "atg", "ref", "annot_tss", "Consensus_5end__TSS", "lr_TSS"]].rename(
    columns={"ref": "refined_TSS", "annot_tss": "v71_TSS", "Consensus_5end__TSS": "experimental_TSS"}
).to_csv(f"{OUT}/Figure5_AD_trimmed_genes.tsv", sep="\t", index=False)
ANC = {"refined": trim.ref.values, "v71": trim.annot_tss.values, "experimental": trim.Consensus_5end__TSS.values}


def meta(sig, anchors, norm):
    out = []
    for c, x, s in zip(trim.chrom.values, anchors, trim.sg.values):
        if not np.isfinite(x): continue
        x = int(x)
        a, b_ = (x - W1, x + W2) if s == 1 else (x - W2 + 1, x + W1 + 1)
        if a < 0 or b_ > CH[c]: continue
        v = sig[c][a:b_]; v = v if s == 1 else v[::-1]
        out.append(v / norm)
    return np.array(out)


pos = np.arange(-W1, W2)
for letter, key, sig, nrm in [("A", "atac", ATAC, GM_ATAC), ("B", "apiap2", AP2, GM_AP2), ("C", "mnase", MN, MED)]:
    df = pd.DataFrame({"bp_from_TSS": pos}); ns = {}
    for an, v in ANC.items():
        M = meta(sig, v, nrm); df[an] = np.nanmean(M, 0); ns[an] = M.shape[0]
    df.to_csv(f"{OUT}/Figure5_{letter}_{key}_metagene.tsv", sep="\t", index=False)
    pd.Series(ns, name="n_genes").to_csv(f"{OUT}/Figure5_{letter}_{key}_n.tsv", sep="\t")

gm, gmr = re.compile("(?=GTGCA)"), re.compile("(?=TGCAC)")
df = pd.DataFrame({"bp_from_TSS": pos}); ns = {}
for an, v in ANC.items():
    out = []
    for c, x, s in zip(trim.chrom.values, v, trim.sg.values):
        if not np.isfinite(x): continue
        x = int(x); sq = SEQ[c][x - W1:x + W2] if s == 1 else rc(SEQ[c][x - W2 + 1:x + W1 + 1])
        z = np.zeros(W1 + W2)
        for m in gm.finditer(sq): z[min(m.start() + 2, W1 + W2 - 1)] += 1
        for m in gmr.finditer(sq): z[min(m.start() + 2, W1 + W2 - 1)] += 1
        out.append(z)
    out = np.array(out); df[an] = out.mean(0) * 1000; ns[an] = out.shape[0]
df.to_csv(f"{OUT}/Figure5_D_motif_GTGCA_per_kb.tsv", sep="\t", index=False)
pd.Series(ns, name="n_genes").to_csv(f"{OUT}/Figure5_D_motif_n.tsv", sep="\t")

# ---------------------------------------------------------------- panel E
bp = pd.read_csv(INPUTS["pairs"], sep="\t")
dft = TR.set_index("gene_id"); H = ["h0", "h8", "h16", "h24", "h32", "h40", "h48"]
L = np.log2(dft[H] + 1)
gb = G.set_index("gene_id")
rng = np.random.default_rng(0); ids = L.index[L.var(axis=1) > 0].values
a = rng.choice(ids, 20000); b_ = rng.choice(ids, 20000)
keep = gb.loc[a, "chrom"].values != gb.loc[b_, "chrom"].values
X = L.loc[a[keep]].values; Y = L.loc[b_[keep]].values
X = X - X.mean(1, keepdims=True); Y = Y - Y.mean(1, keepdims=True)
rnull = (X * Y).sum(1) / np.sqrt((X ** 2).sum(1) * (Y ** 2).sum(1))
E = pd.concat([bp[["orientation", "gene_left", "gene_right", "pearson_DAFT_log2"]].dropna(subset=["pearson_DAFT_log2"])
               .rename(columns={"orientation": "group", "pearson_DAFT_log2": "pearson_r"}),
               pd.DataFrame({"group": "random_interchromosomal", "gene_left": a[keep], "gene_right": b_[keep], "pearson_r": rnull})])
E.to_csv(f"{OUT}/Figure5_E_neighbour_coexpression.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- panels F-G
T = pd.read_csv(INPUTS["pairs_mnase"], sep="\t")
dv = bp[bp.orientation == "divergent"].dropna(subset=["pearson_DAFT_log2"])
P = T.merge(dv[["gene_left", "gene_right", "TSS_left_v71", "TSS_right_v71"]], on=["gene_left", "gene_right"], how="left")
P = P[P.TSS_left_final.notna() & P.TSS_right_final.notna()].copy()
P["dist_v71"] = P.TSS_right_v71 - P.TSS_left_v71
Pm = P[P.dist_v71.notna()].sort_values("TSS_dist_final").reset_index(drop=True)
M = np.full((len(Pm), 2 * WH), np.nan, np.float32)
for i, r in enumerate(Pm.itertuples()):
    mid = int((r.TSS_left_final + r.TSS_right_final) // 2); lo, hi = mid - WH, mid + WH
    if lo < 0 or hi > CH[r.chrom]: continue
    M[i] = np.log2(MN[r.chrom][lo:hi] / MED + 0.05)
np.savez_compressed(f"{OUT}/Figure5_F_mnase_heatmap_matrix.npz", log2_occ_over_median=M.astype(np.float16),
                    bp_from_midpoint=np.arange(-WH, WH))
Pm[["gene_left", "gene_right", "chrom", "TSS_left_final", "TSS_right_final", "TSS_dist_final"]].rename(
    columns={"TSS_left_final": "refined_TSS_left", "TSS_right_final": "refined_TSS_right", "TSS_dist_final": "TSS_distance_bp"}
).to_csv(f"{OUT}/Figure5_F_pair_order.tsv", sep="\t", index=False)

bins = [-1e9, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 1e9]
lab = ["\u2264200", "200\u2013400", "400\u2013600", "600\u2013800", "0.8\u20131k", "1\u20131.5k", "1.5\u20132k", "2\u20133k", ">3k"]
P["distance_bin"] = pd.cut(P.TSS_dist_final, bins, labels=lab)
P["ndr_class"] = np.where(P.n_NDR_interTSS <= 1, "1 NDR", np.where(P.n_NDR_interTSS >= 2, "\u22652 NDRs", "?"))
Gt = P[["gene_left", "gene_right", "TSS_dist_final", "distance_bin", "n_NDR_interTSS", "ndr_class", "r"]].rename(
    columns={"TSS_dist_final": "TSS_distance_bp", "r": "pearson_r"})
Gt.to_csv(f"{OUT}/Figure5_G_pairs_by_distance_and_NDR.tsv", sep="\t", index=False)
s = P[(P.TSS_dist_final >= 300) & (P.TSS_dist_final <= 800) & P.r.notna()]
pG = mannwhitneyu(s.r[s.ndr_class == "1 NDR"], s.r[s.ndr_class == "\u22652 NDRs"]).pvalue

# ---------------------------------------------------------------- panels H-J
DP = pd.read_csv(INPUTS["decomp"], sep="\t"); DC = pd.read_csv(INPUTS["decomp_ctrl"], sep="\t")
DPf = DP[DP.contamination_flag.fillna("") != "exclude (MSP4 read-through)"]
DPf[["pair_id", "gene_id", "desc", "case", "CAGE_dprop", "DAFT_log2FC"]].to_csv(f"{OUT}/Figure5_H_gradeAB_pairs.tsv", sep="\t", index=False)
DC[["pair_id", "gene_id", "case", "CAGE_dprop", "DAFT_log2FC"]].to_csv(f"{OUT}/Figure5_H_singleTSS_controls.tsv", sep="\t", index=False)
foc = {"PF3D7_1311800": "PfA-M1", "PF3D7_1324900": "LDH", "PF3D7_0716800": "eIF3i"}
cols = ["pair_id", "gene_id", "case"] + [f"{d}_prop_{s}" for d in ["CAGE", "U5"] for s in ["ring", "troph", "schizont"]] + \
       [f"DAFT_CPM_{s}" for s in ["ring", "troph", "schizont"]] + ["distal_utr5_len", "proximal_utr5_len"]
I = DPf[DPf.gene_id.isin(foc)][cols].copy(); I.insert(2, "gene_name", I.gene_id.map(foc))
I["order"] = I.gene_id.map({g: i for i, g in enumerate(foc)}); I.sort_values("order").drop(columns="order").to_csv(
    f"{OUT}/Figure5_I_focal_genes.tsv", sep="\t", index=False)
c23 = DPf[DPf.case.str[0].isin(["2", "3"])]
c23[["pair_id", "gene_id", "case", "distal_n_uAUG", "proximal_n_uAUG"]].to_csv(f"{OUT}/Figure5_J_uAUG_pairs.tsv", sep="\t", index=False)
pJ = wilcoxon(c23.distal_n_uAUG, c23.proximal_n_uAUG).pvalue

pd.DataFrame([("G", "Mann-Whitney U, r for 1 vs >=2 NDRs, pairs 300-800 bp apart", len(s), pG),
              ("J", "Wilcoxon signed-rank, distal vs proximal uAUG count", len(c23), pJ)],
             columns=["panel", "test", "n", "p"]).to_csv(f"{OUT}/Figure5_statistics.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- panel geometry (inches) from the draft layout
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize=(7.2, 7.6))
gs = GridSpec(3, 1, figure=fig, height_ratios=[1, 1.05, 1.05], hspace=0.72, left=0.085, right=0.985, top=0.96, bottom=0.1)
r1 = gs[0].subgridspec(1, 4, wspace=0.45); r3 = gs[1].subgridspec(1, 3, wspace=0.62, width_ratios=[0.9, 1, 1.15])
r4 = gs[2].subgridspec(1, 3, wspace=0.55, width_ratios=[1.25, 1.7, 0.75]); sub4 = r4[1].subgridspec(1, 3, wspace=0.12)
spec = [("A", r1[0]), ("B", r1[1]), ("C", r1[2]), ("D", r1[3]), ("E", r3[0]), ("F", r3[1]), ("G", r3[2]),
        ("H", r4[0]), ("I", sub4[0]), ("I", sub4[1]), ("I", sub4[2]), ("J", r4[2])]
rows = []
for k, (pnl, ss) in enumerate(spec):
    bb = ss.get_position(fig)
    rows.append(dict(panel=pnl, order=k, x0_in=bb.x0 * 7.2, y0_in=bb.y0 * 7.6, w_in=bb.width * 7.2, h_in=bb.height * 7.6))
pd.DataFrame(rows).round(4).to_csv(f"{OUT}/Figure5_panel_layout.tsv", sep="\t", index=False)
print("prep done", len(trim), len(Pm), len(DPf), len(c23), round(pG, 3), pJ)
