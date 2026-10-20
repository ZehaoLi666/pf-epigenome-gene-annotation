"""Pool junctions (>=3 unique reads, overhang>=8, canonical motif) and build TPM / tau from featureCounts."""
from pfconfig import HPCC_BIGDATA  # paths: config/paths.yaml
import glob, os, collections, math
import pandas as pd, numpy as np
B = f"{HPCC_BIGDATA}/GSE186109_longread/shortread"
J = B + "/junctions"
STAGES = ["Ring", "Trophozoite", "Schizont", "Gametocyte", "ET", "LT", "Gii", "GV", "Ookinete"]
cols = ["chrom","start","end","strand","motif","annot","uniq","multi","overhang"]
pool = None; per = {}
for s in STAGES:
    d = pd.read_csv(f"{J}/{s}.SJ.out.tab", sep="\t", names=cols)
    per[s] = d.set_index(["chrom","start","end"])["uniq"]
    d = d.set_index(["chrom","start","end","strand","motif","annot"])
    pool = d[["uniq","multi","overhang"]] if pool is None else pool.add(d[["uniq","multi","overhang"]].rename(columns={}), fill_value=0)
pool = pool.reset_index()
ovh = pd.concat([pd.read_csv(f"{J}/{s}.SJ.out.tab", sep="\t", names=cols) for s in STAGES]).groupby(["chrom","start","end"]).overhang.max()
pool = pool.drop(columns="overhang").merge(ovh.reset_index(), on=["chrom","start","end"])
for c in ["uniq","multi","overhang"]: pool[c] = pool[c].astype(int)
pool[cols].to_csv(f"{J}/pooled_all.SJ.out.tab", sep="\t", header=False, index=False)
flt = pool[(pool.uniq >= 3) & (pool.overhang >= 8) & (pool.motif > 0)]
flt[cols].to_csv(f"{J}/pooled_filtered.SJ.out.tab", sep="\t", header=False, index=False)
M = pd.DataFrame(per).fillna(0).astype(int)
M.to_csv(f"{J}/junction_uniq_counts_by_stage.tsv", sep="\t")
print("pooled junctions", len(pool), "filtered", len(flt))
# ---- TPM ----
fc = []
for f in ["pe", "se"]:
    d = pd.read_csv(f"{B}/counts/featurecounts_{f}.txt", sep="\t", comment="#")
    d = d.set_index("Geneid")
    L = d["Length"]
    d = d.iloc[:, 5:]
    d.columns = [os.path.basename(c).replace("_sorted.bam", "") for c in d.columns]
    fc.append(d)
C = pd.concat(fc, axis=1)[STAGES]
C.to_csv(f"{B}/counts/locus_counts_9stages.tsv", sep="\t")
rpk = C.div(L / 1000, axis=0)
TPM = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
TPM["length"] = L
x = np.log2(TPM[STAGES] + 1)
mx = x.max(axis=1)
tau = ((1 - x.div(mx.replace(0, np.nan), axis=0)).sum(axis=1) / (len(STAGES) - 1)).fillna(0)
TPM["tau"] = tau.round(3)
TPM["peak_stage"] = TPM[STAGES].idxmax(axis=1)
TPM["max_TPM"] = TPM[STAGES].max(axis=1).round(2)
TPM["fold_peak_over_median"] = (TPM["max_TPM"] + 1) / (TPM[STAGES].median(axis=1) + 1)
TPM["expressed_any_TPM>=5"] = (TPM[STAGES] >= 5).any(axis=1).astype(int)
TPM.round(3).to_csv(f"{B}/counts/locus_TPM_9stages.tsv", sep="\t")
print("loci", len(TPM), "expressed", int(TPM["expressed_any_TPM>=5"].sum()), "tau>0.8", int((TPM.tau > 0.8).sum()))
assign = []
for f in ["pe", "se"]:
    a = pd.read_csv(f"{B}/counts/featurecounts_{f}.txt.summary", sep="\t", index_col=0)
    a.columns = [os.path.basename(c).replace("_sorted.bam", "") for c in a.columns]
    assign.append(a)
A = pd.concat(assign, axis=1)[STAGES]
A.to_csv(f"{B}/counts/featurecounts_assignment_summary.tsv", sep="\t")
print((A.loc["Assigned"] / A.sum()).round(3).to_dict())
