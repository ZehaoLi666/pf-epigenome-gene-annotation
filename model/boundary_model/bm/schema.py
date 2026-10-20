"""Reference contract: the columns the boundary reference must provide.

One row per BOUNDARY (not per gene). TSS and TTS are separate rows; alternative
ends of the same gene are separate rows sharing gene_id.
"""
from __future__ import annotations
import pandas as pd

REQUIRED = {
    "boundary_id": str,   # unique
    "chrom": str,
    "pos": int,           # 0-based representative coordinate (cluster mode)
    "strand": str,        # '+' or '-'
    "end_type": str,      # 'TSS' or 'TTS'
    "gene_id": str,       # gene / transcription-unit id (lncRNA units allowed)
    "tier": str,          # 'gold' | 'silver' | 'bronze'
}
RECOMMENDED = {
    "iqr_bp": float,          # spread of read ends -> soft-label sigma, sample weight
    "n_reads": int,
    "n_datasets": int,
    "datasets": str,          # comma-separated source ids used to DEFINE the boundary
    "stage_ring": float, "stage_troph": float, "stage_schizont": float, "stage_gam": float,
    "cap_support": float,     # 5' only: cap/TSO/nuclear-total/decay score (higher = genuine)
    "internal_priming": bool, # 3' only
    "linked_cds": str,        # CDS gene id reached by reads from this end, or NA
    "situation": str,         # transcript-graph motif code from the read-end atlas
    "expr_quantile": float,   # 0-1, for stratified evaluation / reweighting
}
TIERS = {"gold", "silver", "bronze"}


def validate_reference(df: pd.DataFrame, chrom_sizes: dict[str, int],
                       validation_sources: set[str] = frozenset()) -> list[str]:
    """Return a list of problems (empty list = valid).

    validation_sources: dataset ids reserved for independent validation
    (e.g. CAGE, poly(A)-anchored 3' ends). They must not appear in `datasets`.
    """
    p: list[str] = []
    miss = [c for c in REQUIRED if c not in df.columns]
    if miss:
        return [f"missing required columns: {miss}"]
    if df["boundary_id"].duplicated().any():
        p.append(f"{df['boundary_id'].duplicated().sum()} duplicated boundary_id")
    bad = ~df["strand"].isin(["+", "-"])
    if bad.any(): p.append(f"{bad.sum()} rows with strand not in +/-")
    bad = ~df["end_type"].isin(["TSS", "TTS"])
    if bad.any(): p.append(f"{bad.sum()} rows with end_type not in TSS/TTS")
    bad = ~df["tier"].isin(TIERS)
    if bad.any(): p.append(f"{bad.sum()} rows with tier not in {sorted(TIERS)}")
    unk = set(df["chrom"]) - set(chrom_sizes)
    if unk: p.append(f"chromosomes not in genome: {sorted(unk)[:5]}")
    sizes = df["chrom"].map(chrom_sizes)
    oob = (df["pos"] < 0) | (df["pos"] >= sizes)
    if oob.any(): p.append(f"{int(oob.sum())} positions outside chromosome bounds")
    if "datasets" in df.columns and validation_sources:
        used = set()
        for s in df["datasets"].dropna():
            used |= {x.strip() for x in str(s).split(",")}
        leak = used & set(validation_sources)
        if leak: p.append(f"LEAKAGE: validation sources used to define labels: {sorted(leak)}")
    if "internal_priming" in df.columns:
        n = int(df.loc[df.end_type == "TTS", "internal_priming"].fillna(False).astype(bool).sum())
        if n: p.append(f"warning: {n} TTS flagged internal_priming still present")
    for c in RECOMMENDED:
        if c not in df.columns:
            p.append(f"note: recommended column absent: {c}")
    return p
