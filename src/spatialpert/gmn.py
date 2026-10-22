"""Geometry-matched null (GMN, AMENDMENT-007).

The null distribution of a knockout's niche statistic comes from pseudo-knockouts made of
control cells with exactly the knockout's lesion geometry (number of lesions and cells per
lesion), placed at random. Assumption: under H0, lesions are seeded independently of position.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.spatial import cKDTree

from spatialpert.effects import between_domain_ss


def niche_ss(resid: np.ndarray, niche: np.ndarray) -> float:
    """Between-niche sum of squares of residuals, cells grouped by their own niche.

    Args:
        resid: Cells x genes residuals.
        niche: Niche per cell.

    Returns:
        T = sum_d n_d ||mean_d - grand mean||^2 summed over genes.
    """
    labs = np.unique(niche)
    sums = np.stack([resid[niche == d].sum(0) for d in labs])
    sizes = np.array([np.sum(niche == d) for d in labs], dtype=np.float64)
    return between_domain_ss(sums, sizes)


def place_geometry(  # noqa: PLR0913 - keyword-only options
    tree: cKDTree,
    ctrl_xy: np.ndarray,
    sizes: Sequence[int],
    rng: np.random.Generator,
    *,
    min_sep: float = 400.0,
    max_tries: int = 5000,
    return_lesions: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Pick control cells forming lesions of the given sizes at random, well-separated seeds.

    Each lesion is a seed plus its nearest unused control cells.

    Args:
        tree: KD-tree on ``ctrl_xy``.
        ctrl_xy: Control cells x 2.
        sizes: Cells per lesion.
        rng: Random generator.
        min_sep: Minimum distance between lesion seeds (µm).
        max_tries: Seed draws before giving up on separation.
        return_lesions: Also return the lesion index of every selected cell.

    Returns:
        Indices into the control cells (concatenated over lesions), and with
        ``return_lesions`` the lesion index per selected cell.
    """
    used = np.zeros(ctrl_xy.shape[0], dtype=bool)
    seeds: list[int] = []
    out: list[np.ndarray] = []
    tries = 0
    for size in sizes:
        while True:
            tries += 1
            s = int(rng.integers(ctrl_xy.shape[0]))
            far = all(np.hypot(*(ctrl_xy[s] - ctrl_xy[t])) > min_sep for t in seeds)
            if (far and not used[s]) or tries > max_tries:
                break
        seeds.append(s)
        _, nn = tree.query(ctrl_xy[s], k=min(int(size) * 3 + 10, ctrl_xy.shape[0]))
        nn = np.atleast_1d(nn)
        nn = nn[~used[nn]][: int(size)]
        used[nn] = True
        out.append(nn)
    idx = np.concatenate(out) if out else np.zeros(0, dtype=int)
    if return_lesions:
        les = (
            np.concatenate([np.full(o.size, j) for j, o in enumerate(out)])
            if out
            else np.zeros(0, dtype=int)
        )
        return idx, les
    return idx


def gmn_pvalue(observed: float, null: np.ndarray) -> float:
    """(1 + #null >= observed) / (1 + M).

    Args:
        observed: Observed statistic.
        null: Null statistics.

    Returns:
        The p-value.
    """
    return float((1 + np.sum(null >= observed)) / (1 + null.size))
