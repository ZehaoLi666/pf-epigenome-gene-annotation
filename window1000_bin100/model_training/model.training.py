#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
End-to-end TSS classifier training:
- Load pos/neg, build labels
- Select TSS-bin features by biology groups
- ANOVA + BH-FDR + effect size filter
- Greedy de-correlation by correlation threshold
- Train XGBoost (CPU default; optional GPU)
- Report metrics and save artifacts

Zehao config defaults baked in; override with CLI if needed.
"""

import os, re, json, argparse, warnings
from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import (
    classification_report, confusion_matrix,
    average_precision_score, roc_auc_score, make_scorer
)
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import f_classif

# --- TEMP PATCH for sklearn 1.6.1 + xgboost 2.1.1 interop ---
import xgboost
def _safe_tags(self):  # minimal tags so sklearn won't choke
    return {"non_deterministic": True}
xgboost.sklearn.XGBClassifier.__sklearn_tags__ = _safe_tags
xgboost.sklearn.XGBRegressor.__sklearn_tags__  = _safe_tags
from xgboost import XGBClassifier


def parse_args():
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[1]

    ap = argparse.ArgumentParser(description="Train XGB on TSS features with ANOVA + de-correlation.")
    ap.add_argument(
        "--pos",
        default=str(repo_root / "window1000_bin100" / "training_features_values" / "positive_values" / "merged_TSS_500bp_20bins.csv"),
    )
    ap.add_argument(
        "--neg",
        default=str(repo_root / "window1000_bin100" / "training_features_values" / "negative_values" / "ATAC_negatives_2k_features_TSS_only.csv"),
    )
    ap.add_argument("--out_dir", default=str(script_dir / "tss_xgb_out"))
    ap.add_argument("--test_size", type=float, default=0.20)
    ap.add_argument("--random_state", type=int, default=42)
    ap.add_argument("--q_fdr_cutoff", type=float, default=0.05)
    ap.add_argument("--eta2_min", type=float, default=0.01)
    ap.add_argument("--corr_thresh", type=float, default=0.95)
    ap.add_argument("--topk_fallback", type=int, default=300)
    ap.add_argument("--n_estimators", type=int, default=2000)
    ap.add_argument("--learning_rate", type=float, default=0.05)
    ap.add_argument("--max_depth", type=int, default=4)
    ap.add_argument("--subsample", type=float, default=0.8)
    ap.add_argument("--colsample_bytree", type=float, default=0.8)
    ap.add_argument("--early_stopping_rounds", type=int, default=200)
    ap.add_argument("--n_splits_cv", type=int, default=5)
    ap.add_argument("--n_jobs", type=int, default=16)
    ap.add_argument("--device", choices=["cpu", "cuda"], default="cpu", help="XGBoost device")
    return ap.parse_args()


# ---------- Helper: PR-AUC scorer for CV ----------
def pr_auc_scorer(est, X, y_true):
    if hasattr(est, "predict_proba"):
        s = est.predict_proba(X)
        s = s[:, 1] if getattr(s, "ndim", 1) > 1 else s
    elif hasattr(est, "decision_function"):
        s = est.decision_function(X)
        s = s[:, 1] if getattr(s, "ndim", 1) > 1 else s
    else:
        s = est.predict(X)
    return average_precision_score(y_true, s)
SCORER_AP = make_scorer(pr_auc_scorer, greater_is_better=True, needs_proba=False)


# ---------- Feature selection: biology groups ----------
GROUP_A = {"H2A.Z", "H2B.Z", "H3K27ac", "H3K9ac", "H3K18ac", "H3K4me3", "H2A.Zac", "ATAC"}   # bins 01–15
GROUP_B = {"H3", "H3K36me3", "H3K4me", "MNase", "H3K27me", "H3K4me2", "H3R17me2", "H3K18me", "H3K4me1"}  # bins 07–15

# e.g. "H3K4me3_TSS_bin07" or "H2A.Zac_TSS_bin01"
PAT = re.compile(r"^(?P<mark>.+?)_TSS_bin(?P<bin>\d+)$")


def keep_column(col: str) -> bool:
    m = PAT.match(col)
    if not m:
        return False
    mark = m.group("mark")
    b = int(m.group("bin"))
    return (mark in GROUP_A and 1 <= b <= 15) or (mark in GROUP_B and 7 <= b <= 15)


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR (safe fallback when statsmodels not available)."""
    m = len(pvalues)
    order = np.argsort(pvalues)
    ranked = np.empty_like(order)
    ranked[order] = np.arange(1, m + 1)
    q = pvalues * m / ranked
    # enforce monotonicity (nonincreasing q-values when ordered by p)
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    q_final = np.empty_like(q_sorted)
    q_final[order] = q_sorted
    return q_final


def anova_filter(X: pd.DataFrame, y: pd.Series, q_fdr_cutoff: float, eta2_min: float, topk_fallback: int):
    imp = SimpleImputer(strategy="median")
    X_imp = pd.DataFrame(imp.fit_transform(X), columns=X.columns, index=X.index)

    F, p = f_classif(X_imp.values, y.values)
    F = np.nan_to_num(F, nan=0.0, posinf=0.0, neginf=0.0)
    p = np.nan_to_num(p, nan=1.0, posinf=1.0, neginf=1.0)

    try:
        from statsmodels.stats.multitest import multipletests
        q = multipletests(p, method="fdr_bh")[1]
    except Exception:
        q = bh_fdr(p)

    k = y.nunique()
    N = len(y)
    eta2 = (F * (k - 1)) / (F * (k - 1) + (N - k) + 1e-12)

    anova_tbl = pd.DataFrame({
        "feature": X.columns,
        "F": F, "p": p, "q": q, "eta2": eta2
    }).sort_values("F", ascending=False).reset_index(drop=True)

    keep = anova_tbl.query("q < @q_fdr_cutoff and eta2 >= @eta2_min")["feature"].tolist()
    if len(keep) == 0:
        keep = anova_tbl.head(min(topk_fallback, len(anova_tbl)))["feature"].tolist()

    return X_imp, anova_tbl, keep


def greedy_decorrelate(X_df: pd.DataFrame, ordered_feats: list, corr_thresh: float) -> list:
    if len(ordered_feats) == 0:
        return []
    corr = X_df[ordered_feats].corr().abs()
    kept = []
    for f in ordered_feats:
        if not kept:
            kept.append(f)
            continue
        if (corr.loc[f, kept].abs() <= corr_thresh).all():
            kept.append(f)
    return kept


def compute_scale_pos_weight(y: pd.Series) -> float:
    n_pos = (y == 1).sum()
    n_neg = (y == 0).sum()
    return float(n_neg / max(n_pos, 1))


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # -------- Load & label
    pos = pd.read_csv(args.pos)
    neg = pd.read_csv(args.neg)
    pos["label"] = 1
    neg["label"] = 0
    data = pd.concat([pos, neg], ignore_index=True)

    # -------- Build X (drop label + optional TSS_coord)
    X_df = data.drop(columns=[c for c in ["label", "TSS_coord"] if c in data.columns], errors="ignore")
    X_df = X_df.select_dtypes(include=["number"])
    print(f"Total numeric features: {X_df.shape[1]}")

    # -------- Select biology-driven TSS-bin columns
    selected_cols = [c for c in X_df.columns if keep_column(c)]
    selected_cols = sorted(selected_cols, key=lambda c: (PAT.match(c).group("mark"), int(PAT.match(c).group("bin"))))
    if len(selected_cols) == 0:
        raise RuntimeError("No columns matched the biology-driven TSS patterns.")

    X_sel = X_df[selected_cols].copy()
    y = data.loc[X_sel.index, "label"].astype(int)
    print(f"Selected {len(selected_cols)}/{X_df.shape[1]} features by biology groups.")
    print(y.value_counts())

    # -------- ANOVA + FDR + effect size
    X_imp, anova_tbl, keep_anova = anova_filter(
        X_sel, y,
        q_fdr_cutoff=args.q_fdr_cutoff,
        eta2_min=args.eta2_min,
        topk_fallback=args.topk_fallback
    )
    ordered_feats = anova_tbl[anova_tbl["feature"].isin(keep_anova)]["feature"].tolist()

    # -------- Greedy de-correlation
    kept = greedy_decorrelate(X_imp, ordered_feats, corr_thresh=args.corr_thresh)
    X_final = X_sel[kept].copy()
    print(f"ANOVA-kept: {len(keep_anova)}  |  After de-correlation: {len(kept)}")

    # -------- Save feature artifacts
    anova_tbl.to_csv(os.path.join(args.out_dir, "anova_table.csv"), index=False)
    pd.Series(keep_anova, name="anova_keep").to_csv(os.path.join(args.out_dir, "features_after_anova.csv"), index=False)
    pd.Series(kept, name="final_features").to_csv(os.path.join(args.out_dir, "features_after_decorrelation.csv"), index=False)

    # -------- Train/val/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X_final, y, test_size=args.test_size, random_state=args.random_state, stratify=y
    )

    # -------- Model config
    spw = compute_scale_pos_weight(y_train)
    print(f"scale_pos_weight={spw:.3f} (neg/pos)")

    model = XGBClassifier(
        objective="binary:logistic",
        n_estimators=args.n_estimators,
        learning_rate=args.learning_rate,
        max_depth=args.max_depth,
        subsample=args.subsample,
        colsample_bytree=args.colsample_bytree,
        reg_lambda=1.0,
        reg_alpha=0.0,
        scale_pos_weight=spw,
        random_state=args.random_state,
        n_jobs=args.n_jobs,
        eval_metric=["aucpr", "auc"],
        tree_method="hist" if args.device == "cpu" else "gpu_hist",
        device=args.device
    )

    # -------- CV on train (PR-AUC)
    skf = StratifiedKFold(n_splits=args.n_splits_cv, shuffle=True, random_state=args.random_state)
    cv_scores = cross_val_score(model, X_train, y_train, cv=skf, scoring=SCORER_AP, n_jobs=1)  # let xgb use n_jobs
    print(f"CV PR-AUC (mean±sd over {args.n_splits_cv} folds): {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
