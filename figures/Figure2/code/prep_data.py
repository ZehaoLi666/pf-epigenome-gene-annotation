"""Gather and derive every input for Figure 2 and write per-panel source tables to ../data.

Inputs (intermediate result tables, read from RESULTS_DIR in config/paths.yaml) used by the original make_fig2.py:

  ml3_end_positives.tsv            ml3_end_positives.tsv  (panel B; was v3tar/work/data/)
  flank_gaps_v3_longread_tiers.tsv flank_gaps_v3_longread_tiers.tsv  (panel B; was v3tar/v3/)
  metagene_profiles.npz            metagene_profiles.npz  (panel C; was fig2v5/)

(make_fig2.py also read fig2v5/input_percentile_final.tsv, but never used it.)

Run:  python prep_data.py   (inputs from RESULTS_DIR), or pass a folder that contains the three files:  python prep_data.py --src DIR
Requires numpy, pandas, scipy, matplotlib.
"""
from pfconfig import result  # paths: config/paths.yaml
import os, sys, shutil, argparse
import numpy as np, pandas as pd
from scipy.stats import trim_mean
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data")
VIDS = {"ml3_end_positives.tsv": "ml3_end_positives.tsv",
        "flank_gaps_v3_longread_tiers.tsv": "flank_gaps_v3_longread_tiers.tsv",
        "metagene_profiles.npz": "metagene_profiles.npz"}


def resolve(src):
    paths = {}
    for fn, vid in VIDS.items():
        if src:
            paths[fn] = os.path.join(src, fn)
        else:
            paths[fn] = result(vid)
        assert os.path.exists(paths[fn]), paths[fn]
    return paths


def layout():
    """Axes rectangles (inches, composite lower-left origin) from the original GridSpec,
    plus the fixed panel canvases (300-dpi pixel grid) and text anchors."""
    W, H = 7.2, 6.6
    fig = plt.figure(figsize=(W, H))
    gs = GridSpec(2, 6, figure=fig, height_ratios=[1.55, 2.25], hspace=0.36, wspace=1.1,
                  left=0.06, right=0.985, top=0.955, bottom=0.06)
    sC = gs[1, :].subgridspec(3, 4, hspace=0.75, wspace=0.32)
    rows = []
    def add(panel, element, kind, x0, y0, w, h):
        rows.append(dict(panel=panel, element=element, kind=kind, order=len(rows),
                         x0_in=x0, y0_in=y0, w_in=w, h_in=h, w_mm=w * 25.4, h_mm=h * 25.4))
    def ax_rect(spec):
        p = spec.get_position(fig); return p.x0 * W, p.y0 * H, p.width * W, p.height * H
    # Canvases: union of each panel's drawn extent in the composite (axes, tick labels,
    # in-panel text, panel C title; panel letters excluded) + ~0.5 mm, snapped to the
    # 300-dpi pixel grid of the composite. Values are pixels from the lower-left corner.
    CANVAS_PX = {"A": (124, 1246, 1404, 1935), "B": (1433, 1198, 2134, 1897), "C": (35, 62, 2160, 1192)}
    for p, (x0, y0, x1, y1) in CANVAS_PX.items():
        add(p, "canvas", "canvas", x0 / 300, y0 / 300, (x1 - x0) / 300, (y1 - y0) / 300)
        if p == "A": add("A", "schematic", "axes", *ax_rect(gs[0, 0:4]))
        if p == "B": add("B", "counts", "axes", *ax_rect(gs[0, 4:6]))
        if p == "C":
            for i in range(12):
                add("C", f"cell{i:02d}", "axes", *ax_rect(sC[i // 4, i % 4]))
            add("C", "title", "text", 0.075 * W, 0.585 * H, 0.0, 0.0)
    plt.close(fig)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--src", default=None); a = ap.parse_args()
    P = resolve(a.src)
    os.makedirs(OUT, exist_ok=True)
    layout().to_csv(os.path.join(OUT, "Figure2_panel_layout.tsv"), sep="\t", index=False)

    # ---------------- panel B: label counts (logic copied from make_fig2.py) -------------
    E = pd.read_csv(P["ml3_end_positives.tsv"], sep="\t")
    FG = pd.read_csv(P["flank_gaps_v3_longread_tiers.tsv"], sep="\t")
    k_ = lambda et, kd, ev=None: int(((E.end_type == et) & (E.kind == kd) &
                                      ((E.end_evidence == ev) if ev is not None else True)).sum())
    rows = [("TSS", "terminal, cap-G", k_("TSS", "terminal", True)),
            ("TSS", "terminal, no cap-G", k_("TSS", "terminal", False)),
            ("TSS", "truncated-model end", k_("TSS", "truncated-model end")),
            ("TSS", "inside CDS (w 0.5)", k_("TSS", "internal")),
            ("TTS", "terminal, poly(A) reads", k_("TTS", "terminal", True)),
            ("TTS", "terminal, no poly(A) reads", k_("TTS", "terminal", False)),
            ("TTS", "truncated-model end", k_("TTS", "truncated-model end")),
            ("TTS", "inside CDS", k_("TTS", "internal")),
            ("Flank", "tier A", int((FG.tier == "A").sum())),
            ("Flank", "tier B", int((FG.tier == "B").sum())),
            ("Flank", "tier C", int((FG.tier == "C").sum())),
            ("Masked gaps", "spanned by reads", int((FG.tier == "spanned").sum()))]
    B = pd.DataFrame(rows, columns=["bar", "segment", "n"])
    B.insert(1, "bar_order", B.bar.map({"TSS": 0, "TTS": 1, "Flank": 2, "Masked gaps": 3}))
    B.insert(3, "segment_order", B.groupby("bar").cumcount())
    B["bar_total"] = B.groupby("bar").n.transform("sum")
    B.to_csv(os.path.join(OUT, "Figure2_B_label_counts.tsv"), sep="\t", index=False)
    # underlying per-label tables (journal source data; not read by panel_B.py)
    shutil.copy(P["ml3_end_positives.tsv"], os.path.join(OUT, "Figure2_B_source_end_positives.tsv"))
    shutil.copy(P["flank_gaps_v3_longread_tiers.tsv"], os.path.join(OUT, "Figure2_B_source_flank_gaps.tsv"))

    # ---------------- panel C: 5 % trimmed-mean metagene profiles ------------------------
    Z = np.load(P["metagene_profiles.npz"])
    ng = {Z[k].shape[0] for k in Z.files}; assert len(ng) == 1
    xb = np.arange(80) + 0.5
    region = np.array(["upstream_1kb"] * 20 + ["gene_body_scaled"] * 40 + ["downstream_1kb"] * 20)
    C = pd.DataFrame({"bin": np.arange(80), "x": xb, "region": region})
    for k in Z.files:
        C[k.replace("_", " ")] = trim_mean(Z[k], 0.05, axis=0)
    C.to_csv(os.path.join(OUT, "Figure2_C_metagene_trimmed_means.tsv"), sep="\t", index=False)
    pd.DataFrame([("n_genes", ng.pop()), ("trim_proportion_each_tail", 0.05),
                  ("bins_upstream_50bp", 20), ("bins_gene_body_scaled", 40), ("bins_downstream_50bp", 20),
                  ("units", "signal / track genome-wide mean")], columns=["key", "value"]
                 ).to_csv(os.path.join(OUT, "Figure2_C_metadata.tsv"), sep="\t", index=False)
    print("wrote", sorted(os.listdir(OUT)))


if __name__ == "__main__":
    main()
