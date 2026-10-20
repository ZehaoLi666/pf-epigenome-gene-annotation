from __future__ import annotations
import numpy as np

# Paper convention: chromosomes 1-5 held out for test. Two further chromosomes
# for validation/early stopping/calibration; never tune on the test set.
def pf_split(chroms: list[str]):
    def num(c):
        s = c.split("_")[1] if "_" in c else c
        s = "".join(ch for ch in s if ch.isdigit())
        return int(s) if s else -1
    test = [c for c in chroms if 1 <= num(c) <= 5]
    val = [c for c in chroms if num(c) in (6, 7)]
    train = [c for c in chroms if c not in test and c not in val and num(c) > 0]
    return train, val, test


def tile(n_bins: int, win: int, step: int):
    starts = list(range(0, max(1, n_bins - win + 1), step))
    if starts[-1] + win < n_bins:
        starts.append(n_bins - win)
    return np.array(starts)
