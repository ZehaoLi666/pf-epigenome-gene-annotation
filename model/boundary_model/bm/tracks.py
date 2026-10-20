"""Feature tracks: manifest -> (n_tracks, n_bins) float32 arrays per chromosome.

Manifest TSV columns: track_id, path, assay, stage, source, role
  role = 'feature' | 'validation'   (validation tracks are never loaded as model input)
"""
from __future__ import annotations
import numpy as np
import pandas as pd

RES = 20  # bp per bin (stage-resolved Pf run)


def load_manifest(path: str) -> pd.DataFrame:
    m = pd.read_csv(path, sep="\t")
    need = {"track_id", "path", "assay", "stage", "source", "role"}
    assert need <= set(m.columns), f"manifest missing {need - set(m.columns)}"
    return m


def feature_tracks(manifest: pd.DataFrame, stage: str | None = None) -> pd.DataFrame:
    m = manifest[manifest.role == "feature"]
    if stage is not None:
        m = m[(m.stage == stage) | (m.stage == "constant")]
    return m.reset_index(drop=True)


def check_partition(manifest: pd.DataFrame, reference_sources: set[str]) -> list[str]:
    """Modality partition guard: no source may be both a feature/label and validation."""
    val = set(manifest.loc[manifest.role == "validation", "source"])
    feat = set(manifest.loc[manifest.role == "feature", "source"])
    probs = []
    if val & feat: probs.append(f"sources both feature and validation: {sorted(val & feat)}")
    if val & reference_sources: probs.append(f"validation sources used in reference: {sorted(val & reference_sources)}")
    return probs


def read_bigwig_binned(path: str, chrom_sizes: dict[str, int], res: int = RES) -> dict[str, np.ndarray]:
    """Mean signal per `res`-bp bin (NaN/no-coverage -> 0)."""
    import pyBigWig
    bw = pyBigWig.open(path)
    have = bw.chroms()
    out = {}
    for c, L in chrom_sizes.items():
        nb = L // res
        if c not in have:
            out[c] = np.zeros(nb, np.float32); continue
        end = min(nb * res, have[c])
        v = np.asarray(bw.values(c, 0, end, numpy=True), dtype=np.float32)
        v = np.pad(v, (0, nb * res - end), constant_values=np.nan)
        # absent intervals = no signal (input-subtracted ChIP clamped at 0; sparse CAGE)
        out[c] = np.nan_to_num(v).reshape(nb, res).mean(1).astype(np.float32)
    bw.close()
    return out


def load_tracks(manifest: pd.DataFrame, chrom_sizes: dict[str, int], res: int = RES) -> dict[str, np.ndarray]:
    per = [read_bigwig_binned(p, chrom_sizes, res) for p in manifest.path]
    return {c: np.stack([t[c] for t in per]).astype(np.float32) for c in chrom_sizes}


def normalise(X: dict[str, np.ndarray], train_chroms: list[str], method: str = "log_z") -> dict[str, np.ndarray]:
    """Per-track normalisation with statistics from TRAINING chromosomes only.

    Apply separately per stage (call once per stage-specific X) so global
    stage-level shifts in a mark do not look like a signature change.
    """
    tr = np.concatenate([X[c] for c in train_chroms], axis=1)
    if method == "log_z":
        tr = np.log1p(np.clip(tr, 0, None))
        mu, sd = tr.mean(1, keepdims=True), tr.std(1, keepdims=True) + 1e-6
        return {c: ((np.log1p(np.clip(v, 0, None)) - mu) / sd).astype(np.float32) for c, v in X.items()}
    if method == "quantile":
        qs = np.linspace(0, 1, 1001)
        ref = np.quantile(tr, qs, axis=1)  # (1001, T)
        out = {}
        for c, v in X.items():
            o = np.empty_like(v)
            for t in range(v.shape[0]):
                o[t] = np.interp(v[t], ref[:, t], qs) * 2 - 1
            out[c] = o
        return out
    raise ValueError(method)
