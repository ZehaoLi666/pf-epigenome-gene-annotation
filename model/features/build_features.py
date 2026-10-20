"""Per-stage feature matrices at 20 bp for the 14 nuclear chromosomes."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import sys, os, json, numpy as np, pandas as pd, pyBigWig
sys.path.insert(0, 'boundary_model')
from bm.tracks import read_bigwig_binned, RES
D = f"{LAB_DIR}/lncRNA project tracks"
STAGES = {"ring": "Ring", "troph": "Trophozoite", "schizont": "Schizont"}
POL = {"ring": "EarlyRing", "troph": "EarlyTroph", "schizont": "LateSchizont"}
bw = pyBigWig.open(f"{D}/histone_marks/H3K4me3_Ring_minusInput.bw")
CS = {c: L for c, L in bw.chroms().items() if c.startswith("Pf3D7_") and c[6:8].isdigit()}
bw.close()
HIST = ["H2AZ", "H2BZ", "H3", "H3K27ac", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3"]
AP2 = sorted({f.split("_")[0] for f in os.listdir(f"{D}/apiap2_tf") if f.endswith(".bw")})
def paths(s):
    S = STAGES[s]; p = {}
    for h in HIST: p[h] = f"{D}/histone_marks/{h}_{S}_minusInput.bw"
    p["ATAC"] = f"{D}/accessibility_nucleosome/ATAC_{S}_minusInput.bw"
    p["MNase"] = f"{D}/accessibility_nucleosome/MNase_{S}_CPM.bw"
    p["HP1"] = f"{D}/machinery/HP1_{S}_minusInput.bw"
    p["PolII_S2P"] = f"{D}/machinery/PolII_{POL[s]}_minusInput.bw"
    for a in AP2: p[a] = f"{D}/apiap2_tf/{a}_{S}_minusInput.bw"
    return p
def rna_paths(s):
    S = STAGES[s]
    pol = {"ring": "Ring_0h", "troph": "Trophozoite_18h", "schizont": "Schizont_36h"}[s]
    return {"GROseq": f"{D}/nascent_transcription/GRO_{S}_mean_normalized.bw",
            "Polysome": f"{D}/translatome/Polysome_{pol}_CPM_final.bw"}
def cage_paths(s):
    return {"CAGE_pos": f"{D}/cage_tss/CTSS_{s if s!='troph' else 'trophozoite'}_pos_CPM.bw",
            "CAGE_neg": f"{D}/cage_tss/CTSS_{s if s!='troph' else 'trophozoite'}_neg_CPM.bw"}
out = {"chrom_sizes": CS, "chrom_tracks": None, "rna_tracks": ["GROseq", "Polysome"]}
for s in STAGES:
    for kind, pf in (("chrom", paths), ("rna", rna_paths), ("cage", cage_paths)):
        P = pf(s)
        for k, v in P.items(): assert os.path.exists(v), v
        per = {k: read_bigwig_binned(v, CS, RES) for k, v in P.items()}
        arr = {c: np.stack([np.abs(per[k][c]) for k in P]).astype(np.float32) for c in CS}
        np.savez(f"data/raw_{kind}_{s}.npz", **arr)
        if kind == "chrom": out["chrom_tracks"] = list(P)
        print(s, kind, len(P), flush=True)
json.dump(out, open("data/meta.json", "w"), indent=1)
